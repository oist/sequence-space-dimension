#!/usr/bin/env bash

# Dimensionality of the sequence space project pipeline

# Project root: /bucket/KondrashovU/seq_space/
# Flash/scratch directory: /flash/KondrashovU/ladaisa/
# Scripts, python binaries (same as src/ in github): /bucket/KondrashovU/seq_space/scripts/
# Scripts, python binaries for simulations (same as src/simulate_evol/ in github): /bucket/KondrashovU/seq_space/simulate_evol/
# Params for simulations folder (same as params/ in github): /bucket/KondrashovU/seq_space/simulate_evol/params/

# Notations: 
# 1 alignment = 1 family = 1 orthogroup (OG)
# dimentionality value for protein family = k coefficient
# AA/aa - amino acids, NT/nt - nucleotides


# ------- NCBI Vertebrates ------

# Alignments (ncbi_vertebrates_data.zip) are provided in the associated Zenodo
# repository - see data/README.md

#----------------1. Vertebrates pairwise distances

#get pairwise distances matrices from ncbi vertebrates alignments
cd alignments_aa/
for i in *fna; do python3 ../../scripts/hamming_distance_multiple_ali_universal.py -a $i \
 -of ../pdistm_matrices/ -d False -aa False ; done

#filter out repeating sequences/points from pdistm matrices (make unique matrices)
head -7292 lists/pdistm_matrices_all_wo_folder.list > lists/pdistm_matrices_all_wo_folder_part1.list
tail -7291 lists/pdistm_matrices_all_wo_folder.list > lists/pdistm_matrices_all_wo_folder_part2.list

sbatch  slurm_scripts/get_pdistm_unique_ncbivert_part1.sh
sbatch slurm_scripts/get_pdistm_unique_ncbivert_part2.sh

#get average pairwise distance
cd /bucket/KondrashovU/seq_space/ncbi_vertebrates
ml python/3.11.4
python ../scripts/get_ld_mean_median_var.py pdistm_matrices_unique/ ncbi_vert_mean_median_var_dist_unique &> get_ld_mean_median_var_pdistm_unique.log &  (wd: /bucket/KondrashovU/seq_space/ncbi_vertebrates)
sed -i 's/_filt_unique//g' allstats_mean_median_var_ld_ncbi_vert_mean_median_var_dist_unique.csv
mv  allstats_mean_median_var_ld_ncbi_vert_mean_median_var_dist_unique.csv tables_sum/

#get max pairwise distance
cd pdistm_matrices_unique
for matr in *pdistm; do awk -F',' '{for(i=1;i<=NF;i++)if($i>max)max=$i} END{sub(/\_filt_unique\.pdistm$/, "", FILENAME); print FILENAME "," max}' ${matr} >> ../vert_max_dist_pdistm_matrices_unique_levenshtein.csv; echo $matr; done &> ../get_max_dist.log &

#get distribution of pairwise distances
cd pdistm_matrices_unique/
sbatch --partition compute -t 10:00:00 -c 1 --mem=100G --output=/flash/KondrashovU/ladaisa/logs/ncbivert/pdistm_hist_0.05.log --job-name=vert_hist005 --wrap "python /bucket/KondrashovU/seq_space/scripts/get_pdist_distribution_plot.py -f /bucket/KondrashovU/seq_space/ncbi_vertebrates/lists/pdistm_matrices_unique_num_of_seq_col_200_seqs_or_more_wo_folder.list -m True -of /flash/KondrashovU/ladaisa/ncbivert/ -fr 0.05"
sbatch --partition compute -t 10:00:00 -c 1 --mem=100G --output=/flash/KondrashovU/ladaisa/logs/ncbivert/pdistm_hist_0.01.log --job-name=vert_hist001 --wrap "python /bucket/KondrashovU/seq_space/scripts/get_pdist_distribution_plot.py -f /bucket/KondrashovU/seq_space/ncbi_vertebrates/lists/pdistm_matrices_unique_num_of_seq_col_200_seqs_or_more_wo_folder.list -m True -of /flash/KondrashovU/ladaisa/ncbivert/ -fr 0.01"


#get distances from pairwise alignemnts for random 5% of the families (616)

cd ncbi_vertebrates/
#get the lis of random families (from the final filtered df - see scripts/analysis/dimension_regression_range_from_data.ipynb)
#move to the taxon folder
mv ../allstats/tmp/random_0.05_og_ids_vert_for_pairali.list .
#rename to correspond to the existing aa ali
sed -i 's/^aa//g' random_0.05_og_ids_vert_for_pairali.list 
sed -i 's/$/_Homo_sapiens_AA_filt_phyid_filtHD.fna.no_duplicates/g' random_0.05_og_ids_vert_for_pairali.list 

#run on the AA sequences from the genomic annotations because MSA are made by codon-aligning the NT sequencing and trimming/removing 
#stop codons and frameshift that depends on the MSA (so just ungapping those sequences wouldn't make them independent)
#get the aa
cd making_ali/sequences_in_all_diam_folder/
mkdir OG_AA_sample
#rename filenames to match aa sequences from annotations
sed 's/_Homo_sapiens_AA_filt_phyid_filtHD.fna.no_duplicates/_Homo_sapiens.faa/g' ../../random_0.05_og_ids_vert_for_pairali.list  > ../../random_0.05_og_ids_vert_for_pairali_orig_aa.list
while read -r f; do unzip -o 18753_vert_OG_AA.zip "$f" -d OG_AA_sample/; done < ../../random_0.05_og_ids_vert_for_pairali_orig_aa.list

#create output dirs
mkdir /flash/KondrashovU/ladaisa/log/ncbivert/pairali/
mkdir /flash/KondrashovU/ladaisa/ncbivert/pairali/
mkdir /flash/KondrashovU/ladaisa/ncbivert/pairali/pdistm_matrices_unique_pairwise_ali
sbatch slurm_scripts/pairwise_ali_unique_pdistm_sample616_ogs.sh 

#move to bucket 
mkdir /bucket/KondrashovU/seq_space/ncbi_vertebrates/pdistm_matrices_unique_pairwise_ali_sample
mkdir /bucket/KondrashovU/seq_space/ncbi_vertebrates/plenm_matrices_unique_pairwise_ali_sample
cd /flash/KondrashovU/ladaisa/ncbivert/pairali/pdistm_matrices_unique_pairwise_ali/
mv *pdistm /bucket/KondrashovU/seq_space/ncbi_vertebrates/pdistm_matrices_unique_pairwise_ali_sample/
mv *plenm /bucket/KondrashovU/seq_space/ncbi_vertebrates/plenm_matrices_unique_pairwise_ali_sample/

#launch dim calculation 
mkdir /flash/KondrashovU/ladaisa/ncbivert/pairali/kcoefs_mult_ali_max_k_range_from_data_pairali
cd /bucket/KondrashovU/seq_space/ncbi_vertebrates/
sed 's/.faa/.pdistm/g' random_0.05_og_ids_vert_for_pairali_orig_aa.list > random_0.05_og_ids_vert_for_pairali_orig_aa_pdistm.list

conda activate prot_evol
sbatch --partition compute -t 10:00:00 -c 1 --mem=20G --output=/flash/KondrashovU/ladaisa/logs/ncbivert/pairali/dimension_reg_range_pairali616.log --job-name=pali_vert_k --wrap "while read -r line; do python3 /bucket/KondrashovU/seq_space/scripts/dimension_steepest_curve_from_matrix.py -f \${line} -i /bucket/KondrashovU/seq_space/ncbi_vertebrates/pdistm_matrices_unique_pairwise_ali_sample/ -of  /flash/KondrashovU/ladaisa/ncbivert/pairali/kcoefs_mult_ali_max_k_range_from_data_pairali/  -w 0.5 -n 20 -r True -p False -s True -m 'max' -l False; done < /bucket/KondrashovU/seq_space/ncbi_vertebrates/random_0.05_og_ids_vert_for_pairali_orig_aa_pdistm.list"

#concatenate
awk 'FNR==1 && NR!=1 {next} {print}' *.coefnr > /bucket/KondrashovU/seq_space/ncbi_vertebrates/vert_pairali_kcoefs_max_k_range_from_data_non_log_0.5_20.coefnr
sed -i 's/_Homo_sapiens//g' /bucket/KondrashovU/seq_space/ncbi_vertebrates/vert_pairali_kcoefs_max_k_range_from_data_non_log_0.5_20.coefnr
sed -i 's/^/aa/g' /bucket/KondrashovU/seq_space/ncbi_vertebrates/vert_pairali_kcoefs_max_k_range_from_data_non_log_0.5_20.coefnr
cd /bucket/KondrashovU/seq_space/
cp ncbi_vertebrates/vert_pairali_kcoefs_max_k_range_from_data_non_log_0.5_20.coefnr allstats/tmp/

#----------------2. Vertebrates dimensionality

#calculate dimension 
#create lists matrices
split -l 1459 --numeric-suffixes lists/pdistm_matrices_all_wo_folder.list lists/pdistm_matrices_all_wo_folder.list
split -l 1459 --numeric-suffixes pdistm_matrices_unique_all_wo_folder.list pdistm_matrices_unique_all_wo_folder.list

#calculate with dimensionality range from the data and non-equal spacing (non-log); unique matrices
sbatch slurm_scripts/dimension_max_k_uniq_range_from_data_0.5_parts10.sh

# win size 0.3 & 0.4
sbatch slurm_scripts/dimension_max_k_uniq_range_from_data_0.3_0.4_parts10.sh

#move to bucket
mv kcoefs_mult_ali_max_k_range_from_dat* /bucket/KondrashovU/seq_space/ncbi_vertebrates/ &

mkdir kcoefs_mult_ali_max_k_range_from_data_unique_0.3/                           
mkdir kcoefs_mult_ali_max_k_range_from_data_unique_0.4/                           
mv /flash/KondrashovU/ladaisa/ncbivert/kcoefs_mult_ali_max_k_range_from_data_unique/*0.3* kcoefs_mult_ali_max_k_range_from_data_unique_0.3/
mv /flash/KondrashovU/ladaisa/ncbivert/kcoefs_mult_ali_max_k_range_from_data_unique/*0.4* kcoefs_mult_ali_max_k_range_from_data_unique_0.4/

#concatenate
cd kcoefs_mult_ali_max_k_range_from_data/
awk 'FNR==1 && NR!=1 {next} {print}' *.coefnr > ../kcoefs_mult_ali_max_k_range_from_data_non_log_0.5_20.coefnr

cd kcoefs_mult_ali_max_k_range_from_data_unique/
awk 'FNR==1 && NR!=1 {next} {print}' *.coefnr > ../kcoefs_mult_ali_max_k_range_from_data_unique_non_log_0.5_20.coefnr

cd kcoefs_mult_ali_max_k_range_from_data_unique_0.3/
awk 'FNR==1 && NR!=1 {next} {print}' *.coefnr > ../kcoefs_mult_ali_max_k_range_from_data_unique_non_log_0.3_20.coefnr
cd kcoefs_mult_ali_max_k_range_from_data_unique_0.4/
awk 'FNR==1 && NR!=1 {next} {print}' *.coefnr > ../kcoefs_mult_ali_max_k_range_from_data_unique_non_log_0.4_20.coefnr

cp kcoefs_mult_ali_max_k_range_from_data_unique_non_log_0.[4,3]_20.coefnr ../allstats/tmp/

#----------------3. Vertebrates dN/dS

#remove stop codons from alignemnts 
sbatch NCBI_gap_stop_codons_for_paml.sh

#calculate dn/ds using two different algorithms in PAML
#get lists of input files 
split -l 18 --numeric-suffixes dnds_NT_ali_list_no_stop_codons_no_paralogs.txt lists_split/dnds_NT_ali_list_no_stop_codons_no_paralogs.txt
cd lists_split/
ls > ../list_of_lists_dnds_NT_ali_list_no_stop_codons_no_paralogs.list

sbatch slurm_scripts/ncbi_vert_dnds_universal_nosplit_ng86_yn00_no_paralogs_14495.sh

#move PAML output files to bucket 
mkdir /bucket/KondrashovU/seq_space/ncbi_vertebrates/dnds/ynout_dir_ng86/
mkdir /bucket/KondrashovU/seq_space/ncbi_vertebrates/dnds/ynout_dir_yn00/
mv /bucket/KondrashovU/seq_space/ncbi_vertebrates/dnds/ynout_dir/ /bucket/KondrashovU/seq_space/ncbi_vertebrates/dnds/ynout_dir_old/
mv ynout_dir/*ng86.ynout /bucket/KondrashovU/seq_space/ncbi_vertebrates/dnds/ynout_dir_ng86/ &
mv ncbi_vert_dnds_* /bucket/KondrashovU/seq_space/ncbi_vertebrates/dnds/
mv ynout_dir/*nostop.ynout /bucket/KondrashovU/seq_space/ncbi_vertebrates/dnds/ynout_dir_yn00/ &

#filter by just removing pairs with ds>=0.8

#ng86
sbatch  ../slurm_scripts/ncbi_vert_14494_parse_ng86_dnds_ynouts.sh
mv /flash/KondrashovU/ladaisa//ncbivert/recalculate_dnds_ng86_yn00/dnds_dn_ds_means_medians_var_fracfilt_ncbi_vert_dnds_ng86_ynout.list.csv /bucket/KondrashovU/seq_space/ncbi_vertebrates/dnds/

#yn00
sbatch  ../slurm_scripts/ncbi_vert_14494_parse_yn00_dnds_ynouts.sh
mv /flash/KondrashovU/ladaisa//ncbivert/recalculate_dnds_ng86_yn00/dnds_dn_ds_means_medians_var_fracfilt_ncbi_vert_dnds_yn00_ynout.list.csv /bucket/KondrashovU/seq_space/ncbi_vertebrates/dnds/

#get table with dn/ds (pdistm are unique)
sed 's/^/aa/g' dnds_dn_ds_means_medians_var_fracfilt_ncbi_vert_dnds_yn00_ynout.list.csv|sed 's/aaOG_id/OG_id/g' > ../tables_sum/dnds_dn_ds_means_medians_var_fracfilt_ncbi_vert_dnds_yn00_ynout_w_aa_ogid.csv

#----------------3. Vertebrates - remove duplicates


#remove duplicates from non-aligned fastas (for now these files aren't used anywhere further)
# removing them from alignemnts is the same as from non-aligned fastas
mkdir /bucket/KondrashovU/seq_space/ncbi_vertebrates/alignments_aa_no_dupl
cd  /bucket/KondrashovU/seq_space/ncbi_vertebrates/alignments_aa
for line in *; do seqkit rmdup -s $line -o ../alignments_aa_no_dupl/$line.no_duplicates; done  &> ../remove_duplicate_vertebrates_14583_fastas.log &

#make all filename suffixes the same 
cd /bucket/KondrashovU/seq_space/ncbi_vertebrates/alignments_aa_no_dupl
mv 12846_Homo_sapiens_AA_muscle_NT_AA_subset_NT_AA_filt_phyid_filtHD.fna.no_duplicates 12846_Homo_sapiens_AA_filt_phyid_filtHD.fna.no_duplicates
mv 16106_Homo_sapiens_AA_muscle_NT_AA_subset_NT_AA_filt_phyid_filtHD.fna.no_duplicates 16106_Homo_sapiens_AA_filt_phyid_filtHD.fna.no_duplicates
mv 4459_Homo_sapiens_filt_muscle_it2_fixed_headers_AA_filt_phyid_filtHD.fna.no_duplicates 4459_Homo_sapiens_AA_filt_phyid_filtHD.fna.no_duplicates
mv 8192_Homo_sapiens_filt_muscle_it2_fixed_headers_AA_filt_phyid_filtHD.fna.no_duplicates 8192_Homo_sapiens_AA_filt_phyid_filtHD.fna.no_duplicates


#----------------4. Vertebrates alignment statistics


#get amino acid usage 

#on unique sequences (identical aa sequences removed to leave only one)
#get the list of aa files with unique sequences
ls *fna > ../lists/ncbi_ver_alignments_aa_no_dupl_14583.list

ml python/3.11.4
sbatch --partition compute -t 24:00:00 -c 1 --mem=20G --output=/flash/KondrashovU/ladaisa/logs/ncbivert/get_vert_usage14583.log --wrap "python /bucket/KondrashovU/seq_space/scripts/get_usage_ali_seq_len_invar_sites_gap_frac.py  -f /bucket/KondrashovU/seq_space/ncbi_vertebrates/lists/ncbi_ver_alignments_aa_no_dupl_14583.list -s '_Homo_sapiens_AA_filt_phyid_filtHD.fna.no_duplicates' -i /bucket/KondrashovU/seq_space/ncbi_vertebrates/alignments_aa_no_dupl/ -o /flash/KondrashovU/ladaisa/ncbivert/"

#move to bucket
mv /flash/KondrashovU/ladaisa/ncbivert/*ali_seq_length_alpha_size_var_invar_sites*.tsv /bucket/KondrashovU/seq_space/ncbi_vertebrates/tables_sum/

#fix og names
cd tables_sum
sed 's/^/aa/g' ncbi_ver_alignments_aa_no_dupl_14583_ali_seq_length_alpha_size_var_invar_sites_keff.tsv |sed 's/aaOG_id/OG_id/g' > ncbi_ver_alignments_aa_no_dupl_14583_ali_seq_length_alpha_size_var_invar_sites_keff_sites_w_aa_ogid.tsv



# ------- Enterobacterales OGs

# Coding sequences from 501 Enterobacterales genomes from NCBI were annotated using PanACoTA pipeline
#the list of genomes with NCBI assembly accession IDs can be found in ./data/metadata/Enterobacterales_genomes501.tsv

#----------------1. Enterobacterales orthologs

#find orthologs

#load modules
ml bioinfo-ugrp-modules  DebianMed/12.0  mmseqs2/14-7e284+ds-1+b2

#create mmseqs database for clustering 
mkdir mmseqs_entrobacterales/
cd mmseqs_entrobacterales/

#select only Enterobacterales species for clustering with mmseqs (finding orthologs)
grep  'Enterobacterales' all_gamma_2801_with_taxonomy.csv|cut -d',' -f1 > mmseqs_entrobacterales/Enterobacterales_gembase_names.list
#get the same list for genes
sed 's/prt$/gen/g' Enterobacterales_gembase_names.list > Enterobacterales_gembase_names_gen.list

#get protein sequences
while read line; do cat Proteins/$line >> mmseqs_entrobacterales/Enterobacterales.prt ; done < mmseqs_entrobacterales/Enterobacterales_gembase_names.list
#get nucleotide sequences
while read line; do cat Genes/$line >> mmseqs_entrobacterales/Enterobacterales.gen ; done < mmseqs_entrobacterales/Enterobacterales_gembase_names_gen.list

#move enterobacteraled database and sequences to the dataset-specific folder
mv /bucket/KondrashovU/seq_space/eggnog/gammaproteo/blast_dnds_orthologs/gammaproteo_annotated_db/mmseqs_entrobacterales/mmseqs_entrobacterales/ /bucket/KondrashovU/seq_space/ncbi_enterobacterales/genomes_db

mmseqs createdb Enterobacterales.prt Enterobacterales.db > createdb.log

#find orthologs
sbatch ../../../slurm_scripts/enterobacterales_mmseqs_find_ogs_id0.5_cov0.8.sh

#get summary
cd /flash/KondrashovU/ladaisa/mmseqs_entrobacterales/
mmseqs createtsv /bucket/KondrashovU/seq_space/eggnog/gammaproteo/blast_dnds_orthologs/gammaproteo_annotated_db/mmseqs_entrobacterales/Enterobacterales.db \
/bucket/KondrashovU/seq_space/eggnog/gammaproteo/blast_dnds_orthologs/gammaproteo_annotated_db/mmseqs_entrobacterales/Enterobacterales.db Enterobacterales_clustered enterobacterales_clustered.tsv

#move to bucket
mkdir /bucket/KondrashovU/seq_space/ncbi_enterobacterales
cd /bucket/KondrashovU/seq_space/ncbi_enterobacterales
mv /flash/KondrashovU/ladaisa/mmseqs_entrobacterales/* .

#get readable result file
mmseqs createseqfiledb /bucket/KondrashovU/seq_space/eggnog/gammaproteo/blast_dnds_orthologs/gammaproteo_annotated_db/mmseqs_entrobacterales/Enterobacterales.db Enterobacterales_clustered enterobacterales_clustered.seq
mmseqs result2flat /bucket/KondrashovU/seq_space/eggnog/gammaproteo/blast_dnds_orthologs/gammaproteo_annotated_db/mmseqs_entrobacterales/Enterobacterales.db /bucket/KondrashovU/seq_space/eggnog/gammaproteo/blast_dnds_orthologs/gammaproteo_annotated_db/mmseqs_entrobacterales/Enterobacterales.db  enterobacterales_clustered.seq enterobacterales_clustered.fasta

#split mmseqs fasta output into separate fastas/orthogroups
mkdir fastas_100seqs/
python ../scripts/split_mmseq.py enterobacterales_clustered.fasta 100 fastas_100seqs/
mkdir fastas_200seqs/
python ../scripts/split_mmseq.py enterobacterales_clustered.fasta 200 fastas_200seqs/

#make list with files 
cd fastas_200seqs/
ls > ../enterobacterales_2605_fastas_200seqs.list
mkdir lists_split
split -l 27 --numeric-suffixes enterobacterales_2605_fastas_200seqs.list lists_split/enterobacterales_2605_fastas_200seqs.list
ls lists_split/enterobacterales_2605_fastas_200seqs.list* > list_of_lists_enterobacterales_2605_fastas_200seqs.list

#move raw mmseqs outputs into a separate directory
mv Enterobacterales_clustered.* mmseqs_raw/
mv enterobacterales_clustered.* mmseqs_raw/


#----------------2. Enterobacterales alignments

#align fastas with muscle
sbatch slurm_scripts/align_aa_enterobacterales200seqs_muscle.sh

#speed up the aligning process by calculating on 470 cores at the same time (instead of 100)

#get the list of already aligned files 
ls -lS /flash/KondrashovU/ladaisa/entero/align_muscle_200seqs/| grep -v 'vuni      0 May' 
nano aligned_fasta_258.list #copy the output (filename column) from the previous command into this file
sed -i 's/\.ali$//g' aligned_fasta_258.list

#get a list of files that haven't been aligned yet 
grep -Fxv -f  aligned_fasta_258.list enterobacterales_2605_fastas_200seqs.list > enterobacterales_2605_fastas_200seqs_remaining_2347.list 

#split into 470 lists
split -l 5 --numeric-suffixes enterobacterales_2605_fastas_200seqs_remaining_2347.list lists_split/enterobacterales_2605_fastas_200seqs_remaining_2347.list 
ls lists_split/enterobacterales_2605_fastas_200seqs_remaining_2347.list* > list_of_lists_enterobacterales_2605_fastas_200seqs_remaining_2347.list 
sbatch slurm_scripts/align_aa_enterobacterales200seqs_muscle_remaining.sh 

#move aligned files from bucket
cd /flash/KondrashovU/ladaisa/entero/align_muscle_200seqs
ls *ali > ../entero_aa_ali_2579.list
while read line; do mv align_muscle_200seqs/$line /bucket/KondrashovU/seq_space/ncbi_enterobacterales/aa_ali_200seqs/; done < /bucket/KondrashovU/seq_space/ncbi_enterobacterales/entero_aa_ali_2579.list &


#----------------3. Enterobacterales pairwise distances

#while all dataset is aligning calculate pairwise distances on a sample 420 families to get an estimate of pairwise distances
#(check that the range of the pairwise distances is comparable to vertebrates on the first 420 families for which alignemnts became available)

#get a list of aligned families (copy the filename column of the output of the following command) and paste it in the file named enterobacterales_sample_420_aa_ali_200seqs.list
ls -lt /flash/KondrashovU/ladaisa/entero/align_muscle_200seqs/| grep -v 'vuni      0 May' 
nano enterobacterales_sample_420_aa_ali_200seqs.list
#calculate levenshtein distances of the sample 420 families
for i in $(seq 1 420); do python3 $INP/scripts/hamming_distance_multiple_ali_universal.py -a $(head -$i $INP/ncbi_enterobacterales/enterobacterales_sample_420_aa_ali_200seqs.list | tail -1) -if $INPF -of $OUTH/ -d True -aa False ; done &
#plot the distribution of pairwise distances in this sample
sed 's/ali/pdistm/g' enterobacterales_sample_420_aa_ali_200seqs.list > enterobacterales_sample_420_pdistm_200seqs.list
cd /flash/KondrashovU/ladaisa/entero/pdistm_matrices_unique_hamming/
#create histogram of pairwise distances (levenshtein)
python /bucket/KondrashovU/seq_space/scripts/get_pdist_distribution_plot.py -f /bucket/KondrashovU/seq_space/ncbi_enterobacterales/enterobacterales_sample_420_pdistm_200seqs.list \
-m True -of /bucket/KondrashovU/seq_space/ncbi_enterobacterales/

#get hist from all
cd /bucket/KondrashovU/seq_space/ncbi_enterobacterales/pdistm_matrices_unique_levenshtein
python /bucket/KondrashovU/seq_space/scripts/get_pdist_distribution_plot.py -f /bucket/KondrashovU/seq_space/ncbi_enterobacterales/enterobacterales_2604_pdistm_200seqs.list \
-m True -of /bucket/KondrashovU/seq_space/ncbi_enterobacterales/ -fr 0.05
python /bucket/KondrashovU/seq_space/scripts/get_pdist_distribution_plot.py -f /bucket/KondrashovU/seq_space/ncbi_enterobacterales/enterobacterales_2604_pdistm_200seqs.list \
-m True -of /bucket/KondrashovU/seq_space/ncbi_enterobacterales/ -fr 0.01

#get list of alignments to calculate pairwise distances
cd aa_ali_200seqs
ls *.ali > ../enterobacterales_2598_aa_ali_200seqs.list
cd ..
split -l 261 --numeric-suffixes enterobacterales_2598_aa_ali_200seqs.list lists_split/enterobacterales_2598_aa_ali_200seqs.list
ls lists_split/enterobacterales_2598_aa_ali_200seqs.list* > list_of_lists_enterobacterales_2598_aa_ali_200seqs.list

#create dirs for outputs
mkdir /flash/KondrashovU/ladaisa/entero/pdistm_matrices_unique_hamming
mkdir /flash/KondrashovU/ladaisa/entero/pdistm_matrices_unique_levenshtein

#get pairwide distance matrices
sbatch  slurm_scripts/entero_2598_200seqs_pairwise_distance_ld_hd.sh 

#the remaining 7 OGs
#get pairwise distances for remaining 7 OGs
while read -r line; do python3 ../scripts/hamming_distance_multiple_ali_universal.py -a $line.ali -i aa_ali_200seqs/ -of pdistm_matrices_unique_levenshtein/ -d True -aa False ; done < enterobacterales_7_remaining_aa_ali_200seqs.list &
while read -r line; do python3 ../scripts/hamming_distance_multiple_ali_universal.py -a $line.ali -i aa_ali_200seqs/ -of pdistm_matrices_unique_hamming/ -d True -aa True -s aa_only ; done < enterobacterales_7_remaining_aa_ali_200seqs.list &

#get average pairwise distance
python ../scripts/get_ld_mean_median_var.py pdistm_matrices_unique_levenshtein/ enterobacterales_unique_levenshtein
sed -i 's/\.fasta//g' allstats_mean_median_var_ld_enterobacterales_unique_levenshtein.csv

#get maximum distance for topological dimension
cd pdistm_matrices_unique_levenshtein
for matr in *pdistm; do awk -F',' '{for(i=1;i<=NF;i++)if($i>max)max=$i} END{sub(/\.fasta\.pdistm$/, "", FILENAME); print FILENAME "," max}' ${matr} >> ../entero_max_dist_pdistm_matrices_unique_levenshtein.csv; echo $matr; done &> ../get_max_dist.log &

#get distances from pairwise alignemnts for random 5% of the families (122)

cd ncbi_enterobacterales/
#get the lis of random families (from the final filtered df - see scripts/analysis/dimension_regression_range_from_data.ipynb)

#move to the taxon folder
mv ../allstats/tmp/random_0.05_og_ids_entero_for_pairali.list .
#rename to correspond to the existing aa ali
sed -i 's/$/.fasta.no_duplicates/g' random_0.05_og_ids_entero_for_pairali.list

#create output dirs
mkdir /flash/KondrashovU/ladaisa/logs/entero/pairali/
mkdir /flash/KondrashovU/ladaisa/entero/pairali/
mkdir /flash/KondrashovU/ladaisa/entero/pairali/pdistm_matrices_unique_pairwise_ali
sbatch slurm_scripts/pairwise_ali_unique_pdistm_sample122_ogs.sh

#move to bucket 
mkdir /bucket/KondrashovU/seq_space/ncbi_enterobacterales/pdistm_matrices_unique_pairwise_ali_sample
mkdir /bucket/KondrashovU/seq_space/ncbi_enterobacterales/plenm_matrices_unique_pairwise_ali_sample
cd /flash/KondrashovU/ladaisa/entero/pairali/pdistm_matrices_unique_pairwise_ali/
mv *pdistm /bucket/KondrashovU/seq_space/ncbi_enterobacterales/pdistm_matrices_unique_pairwise_ali_sample/
mv *plenm /bucket/KondrashovU/seq_space/ncbi_enterobacterales/plenm_matrices_unique_pairwise_ali_sample/

#launch dim calculation 
mkdir /flash/KondrashovU/ladaisa/entero/pairali/kcoefs_mult_ali_max_k_range_from_data_pairali
cd /bucket/KondrashovU/seq_space/ncbi_enterobacterales/
sed 's/.no_duplicates/.pdistm/g' random_0.05_og_ids_entero_for_pairali.list > random_0.05_og_ids_entero_for_pairali_pdistm.list

conda activate prot_evol
sbatch --partition compute -t 10:00:00 -c 1 --mem=20G --output=/flash/KondrashovU/ladaisa/logs/entero/pairali/dimension_reg_range_pairali122.log --job-name=pali_ent_k --wrap "while read -r line; do python3 /bucket/KondrashovU/seq_space/scripts/dimension_steepest_curve_from_matrix.py -f \${line} -i /bucket/KondrashovU/seq_space/ncbi_enterobacterales/pdistm_matrices_unique_pairwise_ali_sample/ -of  /flash/KondrashovU/ladaisa/entero/pairali/kcoefs_mult_ali_max_k_range_from_data_pairali/  -w 0.5 -n 20 -r True -p False -s True -m 'max' -l False; done < /bucket/KondrashovU/seq_space/ncbi_enterobacterales/random_0.05_og_ids_entero_for_pairali_pdistm.list"

#concatenate
awk 'FNR==1 && NR!=1 {next} {print}' *.coefnr > /bucket/KondrashovU/seq_space/ncbi_enterobacterales/entero_pairali_kcoefs_max_k_range_from_data_non_log_0.5_20.coefnr
sed  -i 's/.fasta//g' /bucket/KondrashovU/seq_space/ncbi_enterobacterales/entero_pairali_kcoefs_max_k_range_from_data_non_log_0.5_20.coefnr
cd /bucket/KondrashovU/seq_space/
cp ncbi_enterobacterales/entero_pairali_kcoefs_max_k_range_from_data_non_log_0.5_20.coefnr allstats/tmp/

#----------------4. Enterobacterales dimensionality

#calculate dimensionality 

#calculate dimensionality with regression range from the data and non-equal spacing (non-log); unique matrices
mkdir entero/kcoefs_mult_ali_max_k_range_from_data
mkdir entero/kcoefs_mult_ali_max_k_range_from_data_plots
sbatch ncbi_enterobacterales/slurm_scripts/dimension_max_k_uniq_range_from_data_0.5_ld.sh
mv kcoefs_mult_ali_max_k_range_from_dat* /bucket/KondrashovU/seq_space/ncbi_enterobacterales/ &

#concatenate
cd kcoefs_mult_ali_max_k_range_from_data/
awk 'FNR==1 && NR!=1 {next} {print}' *.coefnr > ../entero_kcoefs_mult_ali_max_k_range_from_data_non_log_0.5_20.coefnr

#window size 0.3 and 0.4 
sbatch slurm_scripts/dimension_max_k_uniq_range_from_data_0.3_0.4_parts10.sh 

#move to bucket 
mv /flash/KondrashovU/ladaisa/entero/kcoefs_mult_ali_max_k_range_from_data/EntOG*0.3* kcoefs_mult_ali_max_k_range_from_data_0.3/ &
mv /flash/KondrashovU/ladaisa/entero/kcoefs_mult_ali_max_k_range_from_data/EntOG*0.4* kcoefs_mult_ali_max_k_range_from_data_0.4/ &

#concatenate
cd kcoefs_mult_ali_max_k_range_from_data_0.3/
awk 'FNR==1 && NR!=1 {next} {print}' *.coefnr > ../entero_kcoefs_mult_ali_max_k_range_from_data_non_log_0.3_20.coefnr

cd kcoefs_mult_ali_max_k_range_from_data_0.4/
awk 'FNR==1 && NR!=1 {next} {print}' *.coefnr > ../entero_kcoefs_mult_ali_max_k_range_from_data_non_log_0.4_20.coefnr

cp entero_kcoefs_mult_ali_max_k_range_from_data_non_log_0.[3,4]_20.coefnr ../allstats/tmp/

#----------------5. Enterobacterales dN/dS

#get nt fasta files
while read -r line; do grep '>' fastas_200seqs/$line > tmp.ids; sed -i 's/^>//g' tmp.ids; seqkit grep -n -f tmp.ids genomes_db/Enterobacterales.gen \
-o nt_fastas_200seqs/$line.nt; done < enterobacterales_2605_fastas_200seqs.list >> get_nt_seqs_for_entogs.log 2>&1 &

#report gaps from aa alignments into nt alignemnts and gap stop codons for PAML
sbatch slurm_scripts/entero_2605_200seqs_report_gaps_aa_to_nt_gap_stops_ali.sh
mv /flash/KondrashovU/ladaisa/entero/NT_alignments_no_stop_codons/ dnds/ &

#report gapt to nt for remaining 7 ali
ml python/3.11.4
for i in $(seq 1 7); do java -jar /apps/unit/KondrashovU/ladaisa/macse_v2.07.jar -prog reportGapsAA2NT \
-align_AA aa_ali_200seqs/$(head -$i enterobacterales_7_remaining_aa_ali_200seqs.list | tail -1).ali -seq nt_fastas_200seqs/$(head -$i enterobacterales_7_remaining_aa_ali_200seqs.list | tail -1).nt \
-out_NT nt_ali_200seqs/$(head -$i enterobacterales_7_remaining_aa_ali_200seqs.list | tail -1).nt.ali; \
python3 /bucket/KondrashovU/seq_space/scripts/gap_stop_codons_for_paml.py -f $(head -$i enterobacterales_7_remaining_aa_ali_200seqs.list | tail -1).nt.ali \
-if nt_ali_200seqs/ -of dnds/NT_alignments_no_stop_codons/ ; done >> remaining_7_report_gaps_stop_codons.log 2>&1

#calculate dn/ds
cd dnds/NT_alignments_no_stop_codons/
ls *ali > ../../entero_2597_nt_ali_wo_stops.list
split -l 4 --numeric-suffixes entero_2597_nt_ali_wo_stops.list lists_split/entero_2597_nt_ali_wo_stops.list
ls entero_2597_nt_ali_wo_stops.list* > ../list_of_lists_entero_2597_nt_ali_wo_stops.list
sbatch /bucket/KondrashovU/seq_space/ncbi_enterobacterales/slurm_scripts/entero_2597_200seqs_dnds_yn00_ng86.sh 

#calculate dn/ds for remaining 7 ali
sbatch slurm_scripts/entero_2597_200seqs_dnds_yn00_ng86_remaining7.sh

#parse dnds 
#NOTE: there are 2 less yn00 parsed files because EntOG0021.fasta.nt_nostop.ynout and EntOG1803.fasta.nt_nostop.ynout had no pairs after filtering

mv /flash/KondrashovU/ladaisa/entero/dnds/ynout_dir/*nostop.ynout /bucket/KondrashovU/seq_space/ncbi_enterobacterales/dnds/ynout_dir_yn00/
mv /flash/KondrashovU/ladaisa/entero/dnds/ynout_dir/*ng86.ynout /bucket/KondrashovU/seq_space/ncbi_enterobacterales/dnds/ynout_dir_ng86/ 
cd ynout_dir_yn00/
ls > ../entero_ynout_yn00_2597.list
cd ../ynout_dir_ng86/
ls > ../entero_ynout_ng86_2597.list
mkdir /flash/KondrashovU/ladaisa/entero/dnds/ynout_dir_yn00_parsed_filtered_seqnames
mkdir /flash/KondrashovU/ladaisa/entero/dnds/ynout_dir_ng86_parsed_filtered_seqnames
sed 's/ynout/ali/g' entero_ynout_yn00_2597.list > entero_ynout_corr_nt_ali_for_seqnames_2597.list
sed -i 's/^/\/bucket\/KondrashovU\/seq_space\/ncbi_enterobacterales\/dnds\/NT_alignments_no_stop_codons\//g' entero_ynout_corr_nt_ali_for_seqnames_2597.list
sbatch ../slurm_scripts/entero_2597_200seqs_parse_dnds_yn00_ng86.sh

#move to bucket 
cd /bucket/KondrashovU/seq_space/ncbi_enterobacterales/dnds
mv /flash/KondrashovU/ladaisa/entero/dnds/dnds_dn_ds_means_medians_var_fracfilt_entero_ynout_ng86_2597.list.csv .
mv /flash/KondrashovU/ladaisa/entero/dnds/dnds_dn_ds_means_medians_var_fracfilt_entero_ynout_yn00_2597.list.csv .

#parse dnds for remaining 7 ali
#prepare lists
cd  /flash/KondrashovU/ladaisa/entero/dnds/ynout_dir/
ls *.nt_nostop.ynout > /bucket/KondrashovU/seq_space/ncbi_enterobacterales/dnds/entero_ynout_yn00_remaining7.list 
ls *.nt_nostop.ng86.ynout > /bucket/KondrashovU/seq_space/ncbi_enterobacterales/dnds/entero_ynout_ng86_remaining7.list
cd /bucket/KondrashovU/seq_space/ncbi_enterobacterales/dnds/
sed -e 's/$/\.nt_nostop\.ali/g' -e 's/^/\/bucket\/KondrashovU\/seq_space\/ncbi_enterobacterales\/dnds\/NT_alignments_no_stop_codons\//g' \
../enterobacterales_7_remaining_aa_ali_200seqs.list > entero_ynout_corr_nt_ali_for_seqnames_remaining7.list

ml python/3.11.4
python3 /bucket/KondrashovU/seq_space/scripts/parse_yn00_output.py \
 -f dnds/entero_ynout_yn00_remaining7.list -id /flash/KondrashovU/ladaisa/entero/dnds/ynout_dir/ -od dnds/ \
 --ynout_dir /flash/KondrashovU/ladaisa/entero/dnds/ynout_dir_yn00_parsed_filtered_seqnames/ -s '.fasta.nt_nostop' -a 'yn00' -l 1 -t True -yo True \
 -sf dnds/entero_ynout_corr_nt_ali_for_seqnames_remaining7.list
 
 python3 /bucket/KondrashovU/seq_space/scripts/parse_yn00_output.py \
 -f dnds/entero_ynout_ng86_remaining7.list -id /flash/KondrashovU/ladaisa/entero/dnds/ynout_dir/ -od dnds/ \
 --ynout_dir /flash/KondrashovU/ladaisa/entero/dnds/ynout_dir_ng86_parsed_filtered_seqnames/ -s '.fasta.nt_nostop.ng86' -a 'ng86' -l 1 -t True -yo True \
 -sf dnds/entero_ynout_corr_nt_ali_for_seqnames_remaining7.list

#move remaining 7 ynouts to bucket
mv /flash/KondrashovU/ladaisa/entero/dnds/ynout_dir/*ng86* ynout_dir_ng86/
mv /flash/KondrashovU/ladaisa/entero/dnds/ynout_dir/*fasta.nt_nostop.ynout ynout_dir_yn00/

#concatenate all yn00 dnds into one final table 
cp dnds_dn_ds_means_medians_var_fracfilt_entero_ynout_yn00_2597.list.csv dnds_dn_ds_means_medians_var_fracfilt_entero_ynout_yn00_2602all.list.csv
tail -7 dnds_dn_ds_means_medians_var_fracfilt_entero_ynout_yn00_remaining7.list.csv >> dnds_dn_ds_means_medians_var_fracfilt_entero_ynout_yn00_2602all.list.csv

cp dnds_dn_ds_means_medians_var_fracfilt_entero_ynout_ng86_2597.list.csv dnds_dn_ds_means_medians_var_fracfilt_entero_ynout_ng86_2604all.list.csv
tail -7 dnds_dn_ds_means_medians_var_fracfilt_entero_ynout_ng86_remaining7.list.csv >> dnds_dn_ds_means_medians_var_fracfilt_entero_ynout_ng86_2604all.list.csv

#move parsed ynouts to bucket
mv /flash/KondrashovU/ladaisa/entero/dnds/ynout_dir_*_parsed_filtered_seqnames/ .

#get list of all 2604 yn00 ynout files
cd /bucket/KondrashovU/seq_space/ncbi_enterobacterales/dnds/ynout_dir_yn00/
ls > ../entero_ynout_yn00_2604.list

#parse yn00 dn/ds while removing all pairs with dn>90 (for cases dn=99 when the program can't calculate it)
ml python/3.11.4
python3 /bucket/KondrashovU/seq_space/scripts/parse_yn00_output.py -f dnds/entero_ynout_yn00_2604.list -id dnds/ynout_dir_yn00/ \ 
-od dnds/ -s '.fasta.nt_nostop' -a 'yn00' -l 1 -t True -dn True &> dnds/parse_yn00_dnds_rem_dn99.log

#parse yn00 dn/ds while removing all pairs with dn>90 (for cases dn=99 when the program can't calculate it)
#and removing all pairs with dn/ds value > 1
ml python/3.11.4
python3 /bucket/KondrashovU/seq_space/scripts/parse_yn00_output.py -f dnds/entero_ynout_yn00_2604.list -id dnds/ynout_dir_yn00/ \
-od dnds/ -s '.fasta.nt_nostop' -a 'yn00' -l 1 -t True -dn True -p True &> dnds/parse_yn00_dnds_rem_dn99_remdnds_more1.log &


#----------------6. Enterobacterales - remove duplicates

#remove duplicates from non-aligned fastas (for now these files aren't used anywhere further)
while read -r line; do seqkit rmdup -s fastas_200seqs/$line -o fastas_200seqs_no_dupl/$line.no_duplicates; done < enterobacterales_2605_fastas_200seqs.list &> remove_duplicate_enterobacterales_2605_fastas_200seqs.log &
#same for nucleotide non-aligned fastas
while read -r line; do seqkit rmdup -s nt_fastas_200seqs/$line -o nt_fastas_200seqs_no_dupl/$line.no_duplicates; done < enterobacterales_2605_nt_fastas_200seqs.list &> remove_duplicate_enterobacterales_2605_nt_fastas_200seqs.log &

#remove duplicates from aligned fastas (amino acids)
cd  /bucket/KondrashovU/seq_space/ncbi_enterobacterales/aa_ali_200seqs
for line in *; do seqkit rmdup -s $line -o ../aa_ali_200seqs_no_dupl/$line.no_duplicates; done  &> ../remove_duplicate_entero_2604_fastas.log &
cd ../aa_ali_200seqs_no_dupl
ls > ../entero_aa_ali_200seqs_no_dupl.list

#get list of duplicated nt sequences in each file to exclude these values from dnds tables
mkdir nt_fastas_200seqs_no_dupl
sed 's/$/\.nt/g' enterobacterales_2605_fastas_200seqs.list > enterobacterales_2605_nt_fastas_200seqs.list
while read -r line; do seqkit rmdup -s nt_fastas_200seqs/$line -o nt_fastas_200seqs_no_dupl/$line.no_duplicates; done < enterobacterales_2605_nt_fastas_200seqs.list &> remove_duplicate_enterobacterales_2605_nt_fastas_200seqs.log &
 
mkdir nt_fastas_200seqs_dupl_seqids/
cd nt_fastas_200seqs
ls > ../enterobacterales_2605_nt_fastas_200seqs.list
 
while read -r line; do grep '>' nt_fastas_200seqs_no_dupl/$line.no_duplicates > tmp_nodup.txt;  grep '>' nt_fastas_200seqs/$line > tmp_dup.txt; grep -Fxv -f  tmp_nodup.txt tmp_dup.txt > nt_fastas_200seqs_dupl_seqids/$line.dupl.ids; done < enterobacterales_2605_nt_fastas_200seqs.list &



#----------------7. Enterobacterales - functional categories + link with other datasets

#link EntOGs with core/cloud + functional categories of COGs
#functional categories are here /bucket/KondrashovU/seq_space/cogs/cog-20.def.tab /bucket/KondrashovU/seq_space/cogs/fun-20.tab
#downloaded from https://ftp.ncbi.nlm.nih.gov/pub/COG/COG2020/data/

#get COG annotations for enterobacterales families from sequences headers 
for file in fastas_200seqs/*; do echo $file >> cog_patterns_aa_fastas_200seqs.txt ; grep -oP '(?<=COG:)COG[0-9]+' $file | sort -u >> cog_patterns_aa_fastas_200seqs.txt; done
sed -i -e 's/\.fasta$//g' -e 's/fastas_200seqs\///g' cog_patterns_aa_fastas_200seqs.txt


#----------------8. Link the Enterobacterales families to the corresponding EggNOG Gammaproteobacteria families

#download and unzip the eggnogs' hmms 
wget http://eggnog5.embl.de/download/eggnog_5.0/per_tax_level/1236/1236_hmms.tar
while read -r line; do tar -xvf 1236_hmms.tar 1236/$line.hmm.gz; done < lists/all_gammaproteo_eggNOGs_aa_w_corresponding_nt_prefixes.list  >> unzip_hmms.log &

#load hmmer 
module load bioinfo-ugrp-modules  DebianMed/12.0 hmmer/3.3.2+dfsg-1

#test on 10 eggnogs 
#get sample
shuf -n 10 lists/all_gammaproteo_eggNOGs_aa_w_corresponding_nt_prefixes.list > lists/all_gammaproteo_eggNOGs_aa_w_corresponding_nt_prefixes_sample10_hmms.list
sed -i 's/$/.hmm.gz/g'  lists/all_gammaproteo_eggNOGs_aa_w_corresponding_nt_prefixes_sample10_hmms.list

#concatenate hmms in one file
while read -r line; do cat 1236_hmm/$line >> hmmer_link_gamma_to_entero/sample10_gamma.hmm; done < lists/all_gammaproteo_eggNOGs_aa_w_corresponding_nt_prefixes_sample10_hmms.list

#run hmmsearch 
hmmsearch --noali --tblout sample10_gamma.tsv -o sample10_gamma.out -E 1e-10 --tformat fasta sample10_gamma.hmm ../../../ncbi_enterobacterales/genomes_db/Enterobacterales.prt &> sample10_gamma_hmmsearch.log

#remove hmm files suffixes to only keep the eggnog family name in query column
sed -i 's/\.faa\.final_tree\.fa//g' sample10_gamma.tsv

#get the table linking each sequence id to EntOG (protein family) - only families with 200 sequences of more
cd /bucket/KondrashovU/seq_space/ncbi_enterobacterales
while read -r line; do grep '>' fastas_200seqs/$line |cut -f1 -d' '|sed -e "s/^>/$line\t/g" -e 's/\.fasta//g' >> entero_200seqs_seqid_to_og.tsv; done < enterobacterales_2605_fastas_200seqs.list &

#run on all dataset
#concatenate hmms in one file
while read -r line; do cat 1236_hmm/$line.hmm.gz >> hmmer_link_gamma_to_entero/all_eggnog_gamma2781.hmm; done < lists/all_gammaproteo_eggNOGs_aa_w_corresponding_nt_prefixes.list &
sed -i 's/\.faa\.final_tree$/\.faa\.final_tree\.fa/g' all_eggnog_gamma2781.hmm

#run hmmer
sbatch /bucket/KondrashovU/seq_space/eggnog/gammaproteo/slurm_scripts/hmmer_link_gamma_to_entero_e-10_all2781.sh

#out of memmory error with 100GB RAM on the first try - only 600 hmms are ready 
#copy these 600 from the output files 
head -1351385 all_eggnog_gamma2781.tsv >all_eggnog_gamma600.tsv
mv all_eggnog_gamma600.tsv /bucket/KondrashovU/seq_space/eggnog/gammaproteo/hmmer_link_gamma_to_entero/

#get the list of remaining hmms 
cd /bucket/KondrashovU/seq_space/eggnog/gammaproteo/
tail -2181 lists/all_gammaproteo_eggNOGs_aa_w_corresponding_nt_prefixes.list > lists/all_gammaproteo_eggNOGs_aa_w_corresponding_nt_prefixes_remaining2181_for_hmms.list

#split 
split -l 546 --numeric-suffixes lists/all_gammaproteo_eggNOGs_aa_w_corresponding_nt_prefixes_remaining2181_for_hmms.list > lists/all_gammaproteo_eggNOGs_aa_w_corresponding_nt_prefixes_remaining2181_for_hmms.list

#get the hmms for each of the four lists
while read -r line; do cat 1236_hmm/$line.hmm.gz >> hmmer_link_gamma_to_entero/all_eggnog_gamma546_remaining1.hmm; done < lists/all_gammaproteo_eggNOGs_aa_w_corresponding_nt_prefixes_remaining2181_for_hmms.list00 &
while read -r line; do cat 1236_hmm/$line.hmm.gz >> hmmer_link_gamma_to_entero/all_eggnog_gamma546_remaining2.hmm; done < lists/all_gammaproteo_eggNOGs_aa_w_corresponding_nt_prefixes_remaining2181_for_hmms.list01 &
while read -r line; do cat 1236_hmm/$line.hmm.gz >> hmmer_link_gamma_to_entero/all_eggnog_gamma546_remaining3.hmm; done < lists/all_gammaproteo_eggNOGs_aa_w_corresponding_nt_prefixes_remaining2181_for_hmms.list02 &
while read -r line; do cat 1236_hmm/$line.hmm.gz >> hmmer_link_gamma_to_entero/all_eggnog_gamma546_remaining4.hmm; done < lists/all_gammaproteo_eggNOGs_aa_w_corresponding_nt_prefixes_remaining2181_for_hmms.list03 &

sed -i 's/\.faa\.final_tree$/\.faa\.final_tree\.fa/g' all_eggnog_gamma546_remaining4.hmm

#run hmmer
sbatch ../slurm_scripts/hmmer_link_gamma_to_entero_e-10_remaining2181_split4.sh

#rename original run
mv all_eggnog_gamma2781.out all_eggnog_gamma2781_first601_orig.out
mv all_eggnog_gamma2781.tsv all_eggnog_gamma2781_first601_orig.tsv

#move everything to bucket 
mv * /bucket/KondrashovU/seq_space/eggnog/gammaproteo/hmmer_link_gamma_to_entero/

sed -i 's/\.faa\.final_tree\.fa//g' all_eggnog_gamma600.tsv 
sed -i 's/\.faa\.final_tree\.fa//g' all_eggnog_gamma546_remaining*.tsv


#----------------9. Count number of sequences 

#unique nucleotides
grep -c '^>' * > ../nt_fastas_no_duplicated_number_of_seq.txt
less nt_fastas_no_duplicated_number_of_seq.txt 
sed 's/\.fasta\.nt\.no_duplicates:/,/g' nt_fastas_no_duplicated_number_of_seq.txt > nt_fastas_no_duplicated_number_of_seq.csv

#unique protein
for afile in pdistm_matrices_unique_levenshtein/*; do grep -c -H '^' $afile ; done  > aa_ali_no_duplicates_num_of_seq.txt  &
sed -e 's/pdistm_matrices_unique_levenshtein\///g' -e 's/\.fasta\.pdistm:/,/g' aa_ali_no_duplicates_num_of_seq.txt > aa_ali_no_duplicates_num_of_seq.csv


#original number of sequences (same for nt and aa)
for file in nt_ali_200seqs/*; do grep -c -H '>' $file ; done  > nt_ali_orig_num_of_seq.txt  &
sed -e 's/nt_ali_200seqs\///g' -e 's/\.fasta\.nt\.ali:/,/g' nt_ali_orig_num_of_seq.txt > nt_aa_ali_orig_num_of_seq.csv


#----------------10. Enterobacterales alignment statistics

cd /bucket/KondrashovU/seq_space/ncbi_enterobacterales/
#on unique sequences (identical aa sequences removed to leave only one)

#get amino acid usage 

#on no duplicates alignments
#get list of 2604 alignemnts
cd aa_ali_200seqs_no_dupl
ls > ../entero_aa_ali_200seqs_no_dupl.list

ml python/3.11.4
sbatch --partition compute -t 3:00:00 -c 1 --mem=20G --output=/flash/KondrashovU/ladaisa/logs/entero/get_entero_usage2604.log --job-name=ent_usage_nd \
--wrap "python /bucket/KondrashovU/seq_space/scripts/get_usage_ali_seq_len_invar_sites_gap_frac.py  -f /bucket/KondrashovU/seq_space/ncbi_enterobacterales/entero_aa_ali_200seqs_no_dupl.list -s '.fasta.ali.no_duplicates' -i /bucket/KondrashovU/seq_space/ncbi_enterobacterales/aa_ali_200seqs_no_dupl/ -o /flash/KondrashovU/ladaisa/entero/"

#on original alignments with duplicated sequences
#get list of 2604 alignemnts
cd aa_ali_200seqs
ls *ali > ../entero_aa_ali_200seqs.list

ml python/3.11.4
sbatch --partition compute -t 3:00:00 -c 1 --mem=20G --output=/flash/KondrashovU/ladaisa/logs/entero/get_entero_usage2604_orig_w_dupl.log --job-name=ent_usage \
--wrap "python /bucket/KondrashovU/seq_space/scripts/get_usage_ali_seq_len_invar_sites_gap_frac.py  -f /bucket/KondrashovU/seq_space/ncbi_enterobacterales/entero_aa_ali_200seqs.list -s '.fasta.ali' -i /bucket/KondrashovU/seq_space/ncbi_enterobacterales/aa_ali_200seqs/ -o /flash/KondrashovU/ladaisa/entero/"

#move to bucket
mv /flash/KondrashovU/ladaisa/entero/*ali_seq_length_alpha_size_var_invar_sites*.tsv /bucket/KondrashovU/seq_space/ncbi_enterobacterales/



# ------- EggNOG v5.0 Gammaproteobacteria ------

#get raw amino acid alignments from EggNOG
wget http://eggnog5.embl.de/download/eggnog_5.0/per_tax_level/1236/1236_raw_algs.tar
#get corresponding nt sequences from NCBI
python get_nt_for_gammaproteo_eggnogs.py

cd /bucket/KondrashovU/seq_space/eggnog/gammaproteo/

#---------------0. Get alignments

#count number of sequences in the nucleotide alignment files
for file in eggNOGs_nt_aligned/*; do grep -c -H '>' $file ; done  > all_gammaproteo_num_of_seq.txt  & 

#count number of sequences in the original amino acid alignment files
for file in 1236_raw_algs/1236/*; do grep -c -H '>' $file ; done  > all_gammaproteo_raw_aa_num_of_seq.txt  &

#turn into csv files
sed 's/\.raw_alg\.faa:/,/g' all_gammaproteo_raw_aa_num_of_seq.txt > all_gammaproteo_raw_aa_num_of_seq.csv

#filter to get only files with 200+ sequences
python3 ../../../scripts/filter_fasta_by_num_of_seq.py allstats_num_of_seq_cut_cogs_4673.txt 200

#same for aa files
python3 ../../scripts/filter_fasta_by_num_of_seq.py all_gammaproteo_raw_aa_num_of_seq.txt 200

#get list of alignment files with 200+ sequences 
tail -1585 all_gammaproteo_eggNOGs_nt_aligned_num_of_seq_200_seqs_or_more.txt |cut -f1 > all_gammaproteo_eggNOGs_nt_aligned_200_seqs_or_more.list

#get the list of aa files 
tail -2784 all_gammaproteo_raw_aa_num_of_seq_200_seqs_or_more.tsv |cut -f1 > all_gammaproteo_eggNOGs_raw_aa_200_seqs_or_more.list

#link gammaproteo eggnogs to cogs 
zgrep 'COG' 1236_annotations.tsv.gz | grep -oP '^.+COG[0-9]+' > eggnog_gamma_cogs_func_link.txt


#get amino acids from nucleotides (non-perfect match between aa and nt sequences creating errors when replacing nt codons in aa alignemnts)
#find eggnog sequences (blast queries in particular) whose nt seqeunce translations don't correspond to their aa sequences
# (as eggnogs don't have their own nt sequences in the database, Olya found the corresponding nucleotides using annotations - see /bucket/KondrashovU/seq_space/eggnog/gammaproteo/get_nt_Olya
# but some nt sequences don't translate exactly to the aa sequences in eggnog which creates a problem when we align sequences and calcualte dnds)

#compare aa sequences from /bucket/KondrashovU/seq_space/eggnog/gammaproteo/blast_dnds_orthologs/aa_with_corresponding_nt/ and nt sequences from 
# /bucket/KondrashovU/seq_space/eggnog/gammaproteo/eggNOGs_nt_aligned/ 

#ungap nt files in  /bucket/KondrashovU/seq_space/eggnog/gammaproteo/eggNOGs_nt_aligned/
#can't use nt files in  /bucket/KondrashovU/seq_space/eggnog/gammaproteo/eggNOGs_nt/ because they have more sequences that files in aa_with_corresponding_nt
sed 's/aa_with_corresponding_nt\///g' blast_dnds_orthologs/aa_with_corresponding_nt_200_seqs_or_more.list > eggNOGs_nt_aligned_aa_with_corresponding_nt_200_seqs_or_more.list
sed -i 's/\.aa_corr_to_nt\.faa/\.fasta/g'  eggNOGs_nt_aligned_aa_with_corresponding_nt_200_seqs_or_more.list
while read -r line; do seqkit seq -g eggNOGs_nt_aligned/$line > eggNOGs_nt_aligned_ungapped/$line; done < eggNOGs_nt_aligned_aa_with_corresponding_nt_200_seqs_or_more.list &

#translate all these nt files
#translate keeping the original first amino acid
while read -r line; do seqkit translate -T 11 eggNOGs_nt_aligned_ungapped/$line > eggNOGs_nt_aligned_ungapped_translated_orig_first_codon/$line; done < eggNOGs_nt_aligned_aa_with_corresponding_nt_200_seqs_or_more.list &
#translate setting the first amino acid to M (methionin) - exactly as the aa sequences from original eggnogs are 
while read -r line; do seqkit translate -T 11 -M eggNOGs_nt_aligned_ungapped/$line > eggNOGs_nt_aligned_ungapped_translated/$line; done < eggNOGs_nt_aligned_aa_with_corresponding_nt_200_seqs_or_more.list &


#find out which nt sequences don't match aa sequences in eggnogs
cd /bucket/KondrashovU/seq_space/eggnog/gammaproteo/
python ../../scripts/compare_nt_translations_to_aa_eggnogs.py blast_dnds_orthologs/aa_with_corresponding_nt_200_seqs_or_more.list eggNOGs_nt_aligned_aa_with_corresponding_nt_200_seqs_or_more.list
#output files are: eggnogs_nt_aa_mismatch_ids.list, eggnogs_nt_aa_mismatch_ids.list 


#realign translated nt 

#get list 
split -l 4 --numeric-suffixes  lists/all_gammaproteo_eggNOGs_nt_translated_1585_and_5remaining_to_align.list lists/all_gammaproteo_eggNOGs_nt_translated_1585_and_5remaining_to_align.list
cd lists/
ls all_gammaproteo_eggNOGs_nt_translated_1585_and_5remaining_to_align.list* > ../list_of_lists_all_gammaproteo_eggNOGs_nt_translated_1585_and_5remaining_to_align.list
mkdir /flash/KondrashovU/ladaisa/eggnoggamma/align_translated_nt_main_1585
sbatch  slurm_scripts/align_translated_nt_over200seqs_1585_and_5rem.sh

#get dnds for remaining 1196 families with less than 200 nt sequences
#create gammaproteo_nt_aligned_ungapped_translated_under200seqs_1196.list list in eggnog_gamma_num_of_seq_and_frac_wo_nucleotides.ipynb
split -l 6 --numeric-suffixes lists/gammaproteo_nt_aligned_ungapped_translated_under200seqs_1196.list lists/gammaproteo_nt_aligned_ungapped_translated_under200seqs_1196.list
cd lists/
ls gammaproteo_nt_aligned_ungapped_translated_under200seqs_1196.list*> ../list_of_lists_gammaproteo_nt_aligned_ungapped_translated_under200seqs_1196.list
mkdir /flash/KondrashovU/ladaisa/eggnoggamma/align_translated_nt
mkdir /flash/KondrashovU/ladaisa/logs/eggnoggamma/align_translated_nt
sbatch  slurm_scripts/align_translated_nt_under200seqs_1196.sh

#move aligned files to bucket
mv /flash/KondrashovU/ladaisa/eggnoggamma/align_translated_nt_main_1585/* eggNOGs_nt_translated_aligned_aa/ &
mv /flash/KondrashovU/ladaisa/eggnoggamma/align_translated_nt/* eggNOGs_nt_translated_aligned_aa/ &

#get list of translated aa from nt (alignments)
cd /bucket/KondrashovU/seq_space/eggnog/gammaproteo/eggNOGs_nt_translated_aligned_aa
ls | sed 's/.ali//g' > ../lists/eggNOGs_nt_translated_aligned_aa_2778.list
split -l 29 --numeric-suffixes  lists/eggNOGs_nt_translated_aligned_aa_2778.list lists/eggNOGs_nt_translated_aligned_aa_2778.list
ls> ../lists/eggNOGs_nt_translated_aligned_aa_2779.list

#original alignments with duplicated sequences
#get list of 2779 alignemnts
cd eggNOGs_nt_translated_aligned_aa
ls  > ../eggNOGs_nt_translated_aligned_aa_2779_orig.list

#----------------1. dN/dS on Gammaproteo



#---1.1 Filter nucleotied alignments for stop codons

#run filtering script 
sbatch  eggnog_gammaprot_gap_stop_codons_for_paml.sh


#---1.2 calculate dN/dS on original EggNOG sequences

#1.2.1 dN/dS on translated nt sequences (latest and main - used for plots)

#report gaps from aa alignments into nt alignemnts and gap stop codons for PAML
mkdir /flash/KondrashovU/ladaisa/eggnoggamma/NT_alignments_no_stop_codons_from_nt/
mkdir /flash/KondrashovU/ladaisa/logs/eggnoggamma/NT_alignments_no_stop_codons_from_nt

sbatch  slurm_scripts/report_gaps_to_nt_from_translated_aligned_2778.sh

mv /flash/KondrashovU/ladaisa/eggnoggamma/NT_alignments_no_stop_codons_from_nt/* /bucket/KondrashovU/seq_space/eggnog/gammaproteo/dnds/NT_alignments_no_stop_codons_from_nt/

cd /bucket/KondrashovU/seq_space/eggnog/gammaproteo/dnds/NT_alignments_no_stop_codons_from_nt/
ls > /bucket/KondrashovU/seq_space/eggnog/gammaproteo/lists/nt_alignments_no_stop_codons_from_nt_2771.list


split -l 10 --numeric-suffixes lists/nt_alignments_no_stop_codons_from_nt_2771.list lists/nt_alignments_no_stop_codons_from_nt_2771.list 

ls lists/nt_alignments_no_stop_codons_from_nt_2771.list* > list_of_lists_nt_alignments_no_stop_codons_from_nt_2771.list

#get remaining 7 OGs
sed 's/.nt_nostop.ali//g' lists/nt_alignments_no_stop_codons_from_nt_2771.list > lists/nt_alignments_no_stop_codons_from_nt_2771_prefixes.list
grep -v -f lists/nt_alignments_no_stop_codons_from_nt_2771_prefixes.list lists/eggNOGs_nt_translated_aligned_aa_2778.list|less
grep -v -f lists/nt_alignments_no_stop_codons_from_nt_2771_prefixes.list lists/eggNOGs_nt_translated_aligned_aa_2778.list > lists/eggNOGs_nt_translated_aligned_aa_remaining_7_to_report_gaps.list
sbatch  slurm_scripts/align_translated_nt_remaining_7.sh 

#calculate dnds
mkdir /flash/KondrashovU/ladaisa/eggnoggamma/yn00_from_nt/ynout_dir
mkdir /flash/KondrashovU/ladaisa/eggnoggamma/yn00_from_nt/logs

sbatch slurm_scripts/eggnog_gamma_from_nt_dnds_yn00_2771.sh
sbatch  slurm_scripts/eggnog_gamma_from_nt_dnds_yn00_remaining_11.sh 

#move to bucket
mkdir  /bucket/KondrashovU/seq_space/eggnog/gammaproteo/dnds/yn00_from_nt/
mv  /flash/KondrashovU/ladaisa/eggnoggamma/yn00_from_nt/ynout_dir/* /bucket/KondrashovU/seq_space/eggnog/gammaproteo/dnds/yn00_from_nt/ &

cd /bucket/KondrashovU/seq_space/eggnog/gammaproteo/dnds/yn00_from_nt/
ls > ../../lists/gammaproteo_from_nt_ynouts_yn00_2775.list

mkdir /bucket/KondrashovU/seq_space/eggnog/gammaproteo/dnds/ynout_dir_yn00_parsed_filtered
mkdir /bucket/KondrashovU/seq_space/eggnog/gammaproteo/dnds/ynout_dir_yn00_parsed_filtered_under1

#parse yn00 outputs 
python3 /bucket/KondrashovU/seq_space/scripts/parse_yn00_output.py -f lists/gammaproteo_from_nt_ynouts_yn00_2775.list -id dnds/yn00_from_nt/ -od dnds/ -s '.fasta.nt_nostop' -a 'yn00' -l 1 -t True -dn True &> dnds/parse_yn00_dnds_2775.log &



#----------------2. Pairwise distances and dimensionality on Gammaproteo

#calculate Levenstein distances on all raw amino acid alignments with >=200 sequences
sbatch gammaproteo_raw_aa2784_pairwise_distance.sh

split -l 1400 --numeric-suffixes lists/egnogg_gamma_raw_pdistm_matrices_all_wo_folder.list lists/egnogg_gamma_raw_pdistm_matrices_all_wo_folder.list

#recalculate dimensionality with regression range from the data and non-equal spacing (non-log); unique matrices
mkdir eggnoggamma/kcoefs_mult_ali_max_k_range_from_data
mkdir eggnoggamma/kcoefs_mult_ali_max_k_range_from_data_plots
#win size 0.5
sbatch eggnog/gammaproteo/slurm_scripts/dimension_max_k_uniq_reg_range_from_data_0.5_parts10.sh
mv  kcoefs_mult_ali_max_k_range_from_dat* /bucket/KondrashovU/seq_space/eggnog/gammaproteo/ &

#win size 0.3 & 0.4
sbatch eggnog/gammaproteo/slurm_scripts/dimension_max_k_uniq_reg_range_from_data_0.3_0.4_parts10.sh

#move to bucket
mv /flash/KondrashovU/ladaisa/eggnoggamma/kcoefs_mult_ali_max_k_range_from_data/*0.3* kcoefs_mult_ali_max_k_range_from_data_0.3/ &
mv /flash/KondrashovU/ladaisa/eggnoggamma/kcoefs_mult_ali_max_k_range_from_data/*0.4* kcoefs_mult_ali_max_k_range_from_data_0.4/ &

#concatenate
cd kcoefs_mult_ali_max_k_range_from_data/
awk 'FNR==1 && NR!=1 {next} {print}' *.coefnr > ../tables_sum/eggnog_gamma_kcoefs_mult_ali_max_k_range_from_data_non_log_0.5_20.coefnr
cd kcoefs_mult_ali_max_k_range_from_data_0.3/
awk 'FNR==1 && NR!=1 {next} {print}' *.coefnr > ../tables_sum/eggnog_gamma_kcoefs_mult_ali_max_k_range_from_data_non_log_0.3_20.coefnr
cd kcoefs_mult_ali_max_k_range_from_data_0.4/
awk 'FNR==1 && NR!=1 {next} {print}' *.coefnr > ../tables_sum/eggnog_gamma_kcoefs_mult_ali_max_k_range_from_data_non_log_0.4_20.coefnr
cp tables_sum/eggnog_gamma_kcoefs_mult_ali_max_k_range_from_data_non_log_0.[4,3]_20.coefnr ../../allstats/tmp/


#get mean, median and var of pairwise distance
python get_ld_mean_median.py all_gammaproteo_pdistm_num_of_seq_200_seqs_or_more_wo_folder.list

#get max distance for topological dimensionality
cd /bucket/KondrashovU/seq_space/eggnog/gammaproteo/pdistm_matrices_unique/
for matr in *pdistm; do awk -F',' '{for(i=1;i<=NF;i++)if($i>max)max=$i} END{sub(/\.raw_alg\.pdistm$/, "", FILENAME); print FILENAME "," max}' ${matr} >> ../gamma_max_dist_pdistm_matrices_unique_levenshtein.csv; echo $matr; done &> ../get_max_dist.log &


#get distribution of pairwise distances
cd pdistm_matrices_unique
python /bucket/KondrashovU/seq_space/scripts/get_pdist_distribution_plot.py -f /bucket/KondrashovU/seq_space/eggnog/gammaproteo/lists/egnogg_gamma_raw_pdistm_matrices_all_wo_folder.list -m True -of /bucket/KondrashovU/seq_space/eggnog/gammaproteo/ -fr 0.05
python /bucket/KondrashovU/seq_space/scripts/get_pdist_distribution_plot.py -f /bucket/KondrashovU/seq_space/eggnog/gammaproteo/lists/egnogg_gamma_raw_pdistm_matrices_all_wo_folder.list -m True -of /bucket/KondrashovU/seq_space/eggnog/gammaproteo/ -fr 0.01


#get distances from pairwise alignemnts for random 5% of the families (134)

cd eggnog/gammaproteo/
#get the lis of random families (from the final filtered df - see scripts/analysis/dimension_regression_range_from_data.ipynb)

#move to the taxon folder
mv ../../allstats/tmp/random_0.05_og_ids_gamma_for_pairali.list .
#rename to correspond to the existing aa ali
sed -i 's/$/.raw_alg.faa/g' random_0.05_og_ids_gamma_for_pairali.list 

#create output dirs
mkdir /flash/KondrashovU/ladaisa/logs/eggnoggamma/pairali/
mkdir /flash/KondrashovU/ladaisa/eggnoggamma/pairali/
mkdir /flash/KondrashovU/ladaisa/eggnoggamma/pairali/pdistm_matrices_unique_pairwise_ali
sbatch  slurm_scripts/pairwise_ali_unique_pdistm_sample134_ogs.sh

#move to bucket 
mkdir /bucket/KondrashovU/seq_space/eggnog/gammaproteo/pdistm_matrices_unique_pairwise_ali_sample
mkdir /bucket/KondrashovU/seq_space/eggnog/gammaproteo/plenm_matrices_unique_pairwise_ali_sample
cd /flash/KondrashovU/ladaisa/eggnoggamma/pairali/pdistm_matrices_unique_pairwise_ali/
mv *pdistm /bucket/KondrashovU/seq_space/eggnog/gammaproteo/pdistm_matrices_unique_pairwise_ali_sample/
mv *plenm /bucket/KondrashovU/seq_space/eggnog/gammaproteo/plenm_matrices_unique_pairwise_ali_sample/
ls *pdistm > random_0.05_og_ids_gamma_for_pairali_pdistm.list
mv random_0.05_og_ids_gamma_for_pairali_pdistm.list /bucket/KondrashovU/seq_space/eggnog/gammaproteo/

#launch dim calculation 
mkdir /flash/KondrashovU/ladaisa/eggnoggamma/pairali/kcoefs_mult_ali_max_k_range_from_data_pairali
cd /bucket/KondrashovU/seq_space/eggnog/gammaproteo/

conda activate prot_evol
sbatch --partition compute -t 10:00:00 -c 1 --mem=20G --output=/flash/KondrashovU/ladaisa/logs/eggnoggamma/pairali/dimension_reg_range_pairali133.log --job-name=pali_gam_k --wrap "while read -r line; do python3 /bucket/KondrashovU/seq_space/scripts/dimension_steepest_curve_from_matrix.py -f \${line} -i /bucket/KondrashovU/seq_space/eggnog/gammaproteo/pdistm_matrices_unique_pairwise_ali_sample/ -of  /flash/KondrashovU/ladaisa/eggnoggamma/pairali/kcoefs_mult_ali_max_k_range_from_data_pairali/  -w 0.5 -n 20 -r True -p False -s True -m 'max' -l False; done < /bucket/KondrashovU/seq_space/eggnog/gammaproteo/random_0.05_og_ids_gamma_for_pairali_pdistm.list"

#concatenate
awk 'FNR==1 && NR!=1 {next} {print}' *.coefnr > /bucket/KondrashovU/seq_space/eggnog/gammaproteo/gamma_pairali_kcoefs_max_k_range_from_data_non_log_0.5_20.coefnr
sed -i 's/.raw_alg//g'  allstats/tmp/gamma_pairali_kcoefs_max_k_range_from_data_non_log_0.5_20.coefnr 
cd /bucket/KondrashovU/seq_space/
cp eggnog/gammaproteo/gamma_pairali_kcoefs_max_k_range_from_data_non_log_0.5_20.coefnr allstats/tmp/

#----------------3. Remove duplicates and count number of unique sequences on Gammaproteobacteria eggnogs

#remove duplicates from aligned fastas (nucleotides) - nt from NCBI
mkdir eggNOGs_nt_aligned_ungapped_unique
while read -r line; do seqkit rmdup -s eggNOGs_nt_aligned_ungapped/$line -o eggNOGs_nt_aligned_ungapped_unique/$line.no_duplicates; done < eggNOGs_nt_aligned_aa_with_corresponding_nt_200_seqs_or_more.list &> remove_duplicate_eggnog_gamma_nt_ungapped_2781_nt_fastas_200seqs.log &

#remove duplicates from aligned fastas (amino acids) - translated nt from NCBI
mkdir eggNOGs_nt_translated_aligned_aa_unique
cd eggNOGs_nt_translated_aligned_aa
for line in *; do seqkit rmdup -s $line -o ../eggNOGs_nt_translated_aligned_aa_unique/$line.no_duplicates; done  &> ../remove_duplicate_eggnog_gamma_nt_ungapped_translated_2779_aa_fastas_200seqs.log &
cd ../eggNOGs_nt_translated_aligned_aa_unique
ls > ../eggNOGs_nt_translated_aligned_aa_2779_unique_no_dupl.list

#count the number of unique nt sequences
for file in eggNOGs_nt_aligned_ungapped_unique/*; do grep -c -H '>' $file ; done  > all_gammaproteo_nt_aligned_ungapped_unique_num_of_seq.txt  & 
sed -e 's/eggNOGs_nt_aligned_ungapped_unique\///g' -e 's/\.fasta\.no_duplicates:/,/g' all_gammaproteo_nt_aligned_ungapped_unique_num_of_seq.txt > all_gammaproteo_nt_aligned_ungapped_unique_num_of_seq.csv

#count number of unique amino acid sequences from pdistm matrices (duplicates filtered during pdistm calculation)
for afile in pdistm_matrices_unique/*; do grep -c -H '^' $afile ; done  > all_gammaproteo_pdistm_num_of_seq.txt  &
sed -e 's/pdistm_matrices_unique\///g' -e 's/\.raw_alg\.pdistm:/,/g' all_gammaproteo_pdistm_num_of_seq.txt > all_gammaproteo_pdistm_num_of_seq.csv

#----------------4. Gammaproteobacteria eggnogs alignment statistics & Evolution parameters

cd /bucket/KondrashovU/seq_space/eggnog/gammaproteo/

#get amino acid usage 
#original EggNOG amino acids
#on no duplicates alignments
ml python/3.11.4
sbatch --partition compute -t 3:00:00 -c 1 --mem=20G --output=/flash/KondrashovU/ladaisa/logs/eggnoggamma/get_gamma_orig_aa_usage2784.log --job-name=gamma_usage_origaa \
--wrap "python /bucket/KondrashovU/seq_space/scripts/get_usage_ali_seq_len_invar_sites_gap_frac.py  -f /bucket/KondrashovU/seq_space/eggnog/gammaproteo/lists/all_gammaproteo_eggNOGs_raw_aa_200_seqs_or_more.list -s '.raw_alg.faa' -i /bucket/KondrashovU/seq_space/eggnog/gammaproteo/1236_raw_algs/1236/ -o /flash/KondrashovU/ladaisa/eggnoggamma/"

#move to bucket 
cd /flash/KondrashovU/ladaisa/eggnoggamma
mv  *invar*tsv /bucket/KondrashovU/seq_space/eggnog/gammaproteo/tables_sum/



# ------- COG prokaryotes - full genes ------

#sequences were downloaded from https://ftp.ncbi.nlm.nih.gov/pub/COG/COG2020/data/
wget https://ftp.ncbi.nlm.nih.gov/pub/COG/COG2020/data/cog-20.fa.gz
gunzip cog-20.fa.gz

#----------------1. COGs filtering

#count number of unique sequences (size of the pdistm_unique matrices) 
for line in sequences/*fa; do grep -c '>'  $line >> cogs_full_num_of_seq.csv; done

cd sequences/
ls *fa > ../lists/cogs_seq_full_wo_folder.list
cd ../
paste lists/cogs_seq_full_wo_folder.list cogs_full_num_of_seq.csv  -d "," > cogs4877_num_of_seq.csv

#filter by size (keep only ones with 200+ sequences)
sed -i 's/,/:/g' cogs4877_num_of_seq.csv
python3 ../../scripts/filter_fasta_by_num_of_seq.py cogs4877_num_of_seq.csv 200
tail -2959 cogs4877_num_of_seq_200_seqs_or_more.tsv| cut -f 1 > lists/cogs2959_num_of_seq_200_seqs_or_more.list


#----------------2. COGs Pairwise alignments and distances  

#split in lists of 20 for pairwise alignemnt 
split -l 30 --numeric-suffixes lists/cogs2959_num_of_seq_200_seqs_or_more.list lists/cogs_num_of_seq_200_seqs_or_more_split.list
cd lists/
ls cogs_num_of_seq_200_seqs_or_more_split.list* > list_of_lists_cogs_num_of_seq_200_seqs_or_more_split.list

#split in lists of 10 for pairwise alignemnt 
split -l 10 --numeric-suffixes lists/cogs2959_num_of_seq_200_seqs_or_more.list lists/cogs_num_of_seq_200_seqs_or_more_split.list
cd lists/
ls cogs_num_of_seq_200_seqs_or_more_split.list* > list_of_lists_cogs_num_of_seq_200_seqs_or_more_split.list

#calculate pairwise distances (run this as many times as needed - it writes partial progress and can resrart from the place where it finished)
sbatch slurm_scripts/pairwise_ali_hamming_distance_cogs_continue.sh

#remove finished from the list and use more cores
cd /flash/KondrashovU/ladaisa/cogs/full/muscle/pdistm_matrices_unique_pairwise_ali/
ls *pdistm |sed 's/.pdistm/.fa/g' >../done07092025.list
grep -vxF -f ../done07092025.list /bucket/KondrashovU/seq_space/cogs/full/lists/cogs2959_num_of_seq_200_seqs_or_more.list > /bucket/KondrashovU/seq_space/cogs/full/lists/cogs2959_num_of_seq_200_seqs_or_more_remaining2399.list

cd lists/
split -l 4 --numeric-suffixes cogs2959_num_of_seq_200_seqs_or_more_remaining2399.list cogs2959_num_of_seq_200_seqs_or_more_remaining2399_split.list
ls cogs2959_num_of_seq_200_seqs_or_more_remaining2399_split.list* > list_of_lists_cogs2959_num_of_seq_200_seqs_or_more_remaining2399_split.list

sbatch /bucket/KondrashovU/seq_space/cogs/full/slurm_scripts/pairwise_ali_hamming_distance_cogs_continue_rem.sh

#repeat 
cd /flash/KondrashovU/ladaisa/cogs/full/muscle/pdistm_matrices_unique_pairwise_ali/
ls *pdistm |sed 's/.pdistm/.fa/g' >../done12092025.list
cd ../
# add manually these 4 families (not enough memory - to not rerun them again with the same small RAM): COG3170.fa, COG3319.fa, COG5276.fa, COG5281.fa
grep -vxF -f done12092025.list /bucket/KondrashovU/seq_space/cogs/full/lists/cogs2959_num_of_seq_200_seqs_or_more.list > /bucket/KondrashovU/seq_space/cogs/full/lists/cogs2959_num_of_seq_200_seqs_or_more_remaining1830.list
cd /bucket/KondrashovU/seq_space/cogs/full/lists/
split -l 3 --numeric-suffixes cogs2959_num_of_seq_200_seqs_or_more_remaining1830.list cogs2959_num_of_seq_200_seqs_or_more_remaining1830_split.list
ls cogs2959_num_of_seq_200_seqs_or_more_remaining1830_split.list* > list_of_lists_cogs2959_num_of_seq_200_seqs_or_more_remaining1830_split.list
sbatch  slurm_scripts/pairwise_ali_hamming_distance_cogs_continue_rem.sh 

#repeat 
cd /flash/KondrashovU/ladaisa/cogs/full/muscle/pdistm_matrices_unique_pairwise_ali/
ls *pdistm |sed 's/.pdistm/.fa/g' >../done22092025.list
cd ../
# add manually these 4 families (not enough memory - to not rerun them again with the same small RAM): COG3170.fa, COG3319.fa, COG5276.fa, COG5281.fa
grep -vxF -f done22092025.list /bucket/KondrashovU/seq_space/cogs/full/lists/cogs2959_num_of_seq_200_seqs_or_more.list > /bucket/KondrashovU/seq_space/cogs/full/lists/cogs2959_num_of_seq_200_seqs_or_more_remaining1343.list
cd /bucket/KondrashovU/seq_space/cogs/full/lists/
split -l 2 --numeric-suffixes cogs2959_num_of_seq_200_seqs_or_more_remaining1343.list cogs2959_num_of_seq_200_seqs_or_more_remaining1343_split.list
ls cogs2959_num_of_seq_200_seqs_or_more_remaining1343_split.list* > list_of_lists_cogs2959_num_of_seq_200_seqs_or_more_remaining1343_split.list
sbatch  slurm_scripts/pairwise_ali_hamming_distance_cogs_continue_rem.sh 

#repeat
cd /flash/KondrashovU/ladaisa/cogs/full/muscle/pdistm_matrices_unique_pairwise_ali/
ls *pdistm |sed 's/.pdistm/.fa/g' >../done30102025.list
cd ../
# add manually these 4 families (not enough memory - to not rerun them again with the same small RAM): COG3170.fa, COG3319.fa, COG5276.fa, COG5281.fa
grep -vxF -f done30102025.list /bucket/KondrashovU/seq_space/cogs/full/lists/cogs2959_num_of_seq_200_seqs_or_more.list > /bucket/KondrashovU/seq_space/cogs/full/lists/cogs2959_num_of_seq_200_seqs_or_more_remaining743.list
cd /bucket/KondrashovU/seq_space/cogs/full/lists/
sbatch  slurm_scripts/pairwise_ali_hamming_distance_cogs_continue_rem.sh 

sbatch slurm_scripts/pairwise_ali_hamming_distance_cogs_continue_rem_rest76.sh
sbatch cogs/full/slurm_scripts/pairwise_ali_hamming_distance_cogs_rm_long_rest.sh
sbatch cogs/full/slurm_scripts/pairwise_ali_hamming_distance_cogs_rm_long_rest530.sh


#get sample of pairwise distances for the distribution plot
sed 's/^/\/flash\/KondrashovU\/ladaisa\/cogs\/full\/muscle\/pdistm_matrices_unique_pairwise_ali\//g' pairali_done12112025.list > pairali_done12112025_wpath.list
sbatch --partition compute -t 10:00:00 -c 1 --mem=100G --output=/flash/KondrashovU/ladaisa/logs/cogs/full/pdistm_hist_0.05.log --job-name=cogsfull_hist005 --wrap "python /bucket/KondrashovU/seq_space/scripts/get_pdist_distribution_plot.py -f /bucket/KondrashovU/seq_space/cogs/full/pairali_done12112025_wpath.list -m True -of /flash/KondrashovU/ladaisa/cogs/full/ -fr 0.05"
mv /flash/KondrashovU/ladaisa/cogs/full/pairali_done12112025_wpath* .
cp pairali_done12112025_wpath_5_pairs_hist.pickle ../../compare_all_datasets/



#get max distance
cd /flash/KondrashovU/ladaisa/cogs/full/muscle/pdistm_matrices_unique_pairwise_ali/
for matr in *pdistm; do awk -F',' '{for(i=1;i<=NF;i++)if($i>max)max=$i} END{sub(/\.pdistm$/, "", FILENAME); print FILENAME "," max}' ${matr} >> /bucket/KondrashovU/seq_space/cogs/full/cogs_full_pairali_max_dist_pdistm_matrices_unique_levenshtein.csv; echo $matr; done &> /bucket/KondrashovU/seq_space/cogs/full/get_max_dist.log &

#get max distance for next 308 COGs
cd /flash/KondrashovU/ladaisa/cogs/full/muscle/pdistm_matrices_unique_pairwise_ali/
while read -r matr; do awk -F',' '{for(i=1;i<=NF;i++)if($i>max)max=$i} END{sub(/\.pdistm$/, "", FILENAME); print FILENAME "," max}' ${matr} >> /bucket/KondrashovU/seq_space/cogs/full/cogs_full_pairali_max_dist_pdistm_matrices_unique_levenshtein_308.csv; echo $matr; done < /bucket/KondrashovU/seq_space/cogs/full/pairali_done12112025_last308.list &> /bucket/KondrashovU/seq_space/cogs/full/get_max_dist_308.log &

# get mean, median and variance of LD for all 2520 COGs
python ../../scripts/get_ld_mean_median_var.py /flash/KondrashovU/ladaisa/cogs/full/muscle/pdistm_matrices_unique_pairwise_ali/ cogs_full_pairali_unique_levenshtein_2520 > get_ld_mean_median_var_2520.log &
sed -i 's/\.fasta//g' allstats_mean_median_var_ld_cogs_full_pairali_unique_levenshtein_2520.csv



#----------------3. COGs dimensionality

#calculate dimensionality from the finished 2212 COGs
ls *pdistm |sed 's/.pdistm/.fa/g' >../done30102025.list
cp /flash/KondrashovU/ladaisa/cogs/full/muscle/done30102025.list pairali_done30102025.list
sed -i 's/.fa/.pdistm/g' pairali_done30102025.list
mkdir /flash/KondrashovU/ladaisa/cogs/full/kcoefs_mult_ali_max_k_range_from_data_pairali
mkdir /flash/KondrashovU/ladaisa/cogs/full/kcoefs_mult_ali_max_k_range_from_data_pairali_plots
sbatch  slurm_scripts/dim_range_from_data_from_pairwise_ali.sh

#concatenate
cd /flash/KondrashovU/ladaisa/cogs/full/kcoefs_mult_ali_max_k_range_from_data_pairali
awk 'FNR==1 && NR!=1 {next} {print}' *.coefnr > /bucket/KondrashovU/seq_space/cogs/full/cogs_full_pairali_kcoefs_mult_ali_max_k_range_from_data_non_log_0.5_20.coefnr

#get dim for the next 308 COGs
cd /flash/KondrashovU/ladaisa/cogs/full/muscle/pdistm_matrices_unique_pairwise_ali/
ls *pdistm > ../done12112025.list
cd /bucket/KondrashovU/seq_space/cogs/full/
cp /flash/KondrashovU/ladaisa/cogs/full/muscle/done12112025.list pairali_done12112025.list
sed -i 's/.fa/.pdistm/g' pairali_done12112025.list 
grep -v -f pairali_done30102025.list pairali_done12112025.list > pairali_done12112025_last308.list
sbatch  slurm_scripts/dim_range_from_data_from_pairwise_ali.sh

#concatenate next 308 COGs
cd /flash/KondrashovU/ladaisa/cogs/full/kcoefs_mult_ali_max_k_range_from_data_pairali
awk 'FNR==1 && NR!=1 {next} {print}' *.coefnr > /bucket/KondrashovU/seq_space/cogs/full/cogs_full_pairali_kcoefs_mult_ali_max_k_range_from_data_non_log_0.5_20_308.coefnr


#----------------4. COGs statistics


#get the mean and median of the sequence length
python ../../scripts/get_mean_median_seq_len_from_fasta.py done30102025_wo_big.list sequences/ cogs_full_pairali_mean_median_seq_len_orig.csv > get_mean_median_seq_len_from_fasta_2212.log &

#get the mean and median of the sequence length (next 308)
python ../../scripts/get_mean_median_seq_len_from_fasta.py pairali_done12112025_last308_fasta.list sequences/ cogs_full_pairali_mean_median_seq_len_orig_308.csv COG1196,COG2931,COG3170,COG3210,COG3319,COG3321,COG5276,COG5281 10000 > get_mean_median_seq_len_from_fasta_308.log &

#get number of unique sequences (rows in pdistm)
cd /flash/KondrashovU/ladaisa/cogs/full/muscle/pdistm_matrices_unique_pairwise_ali/
for f in *.pdistm; do echo "$f,$(wc -l < "$f")"; done > ../num_of_seq_unique_cogs2212.csv
mv  ../num_of_seq_unique_cogs2212.csv /bucket/KondrashovU/seq_space/cogs/full/
cd /bucket/KondrashovU/seq_space/cogs/full/
sed -i 's/.pdistm//g' num_of_seq_unique_cogs2212.csv 

#get number of unique sequences (rows in pdistm) (next 308)
cd /flash/KondrashovU/ladaisa/cogs/full/muscle/pdistm_matrices_unique_pairwise_ali/
while read -r f; do echo "$f,$(wc -l < "$f")"; done < /bucket/KondrashovU/seq_space/cogs/full/pairali_done12112025_last308.list > /bucket/KondrashovU/seq_space/cogs/full/num_of_seq_unique_cogs308.csv

cd /bucket/KondrashovU/seq_space/cogs/full/
sed -i 's/.pdistm//g' num_of_seq_unique_cogs308.csv 


#COGs with too long sequences: COG3170, COG3319, COG5276, COG5281, COG1196, COG2931, COG3210, COG3321
#COGs with too many sequences (20000+): COG1028, COG0583
#COGs with too many sequences (<20000):
COG2197.fa      10241
COG2226.fa      10804
COG2207.fa      12177
COG1595.fa      12319
COG0596.fa      13171
COG2202.fa      13407
COG2814.fa      16066
COG0438.fa      16297
COG0745.fa      17683
COG1309.fa      18838
COG0642.fa      18976
COG1028.fa      20713
COG0583.fa      21908



# -------- Simulating protein families to study the parameters influencing the dimensionality and dn/ds -----------

# params/ - folder with the lists of simulation parameters (inputs for simulate_protein_family*.sh files)


### Simulating branch length dependence on tree depth (+ epistasis)

# Simulations for tree and notree and gamma 0,1,10,100 each AND for F 0.07-0.7):
#here we also record the distance matrix and calculate dim on leaves, ancestor nodes, and leaves with ancestor nodes

#run short/normal leaves and notree with 100GB RAM 
sbatch simulate_protein_family_evol_diff_notree_tree_cont_f_normshortleaves_w_anc_distm.sh
#run longleaves and notree with 450GB RAM 
sbatch simulate_protein_family_evol_diff_notree_tree_cont_f_longleaves_w_anc_distm.sh

#notree with 100GB RAM failed with OOM error, restart them with 450GB
sbatch simulate_evol/simulate_protein_family_evol_diff_notree_tree_cont_f_notree_w_anc_distm.sh

#about 100-160 simulations ran with previous runs, get the remaining ones with 10 simulations per run
sbatch simulate_protein_family_evol_diff_notree_tree_cont_f_notree0.7_10reps_w_anc_distm.sh 


mkdir /bucket/KondrashovU/seq_space/simulate_evol/data/uniform_f_w_dist
mkdir /bucket/KondrashovU/seq_space/simulate_evol/data/uniform_f_w_dist/csv
mkdir /bucket/KondrashovU/seq_space/simulate_evol/data/uniform_f_w_dist/fasta
mkdir /bucket/KondrashovU/seq_space/simulate_evol/data/uniform_f_w_dist/trees
mkdir /bucket/KondrashovU/seq_space/simulate_evol/data/uniform_f_w_dist/logs
mkdir /bucket/KondrashovU/seq_space/simulate_evol/data/uniform_f_w_dist/pdistm

cd /flash/KondrashovU/ladaisa/simulate_evol/uniform_frac_w_dist

#get list of all fastas
ls > ../final_all_fastas_5244.list
cd uniform_f_w_dist/
mkdir lists
split -l 200 --numeric-suffixes final_all_fastas_5244.list lists/final_all_fastas_5244.list
cd lists/
ls > ../list_of_list_final_all_fastas_5244.list

mv *csv /bucket/KondrashovU/seq_space/simulate_evol/data/uniform_f_w_dist/csv/

#remove duplicates from pdistm matrices, rename them for Betti0 curves calculation and zip 
sbatch filter_duplicated_seqs_from_simu_pdistm_rename_and_zip.sh

#get lists of trees for each condition
ls *notree_1_gammaFitAll*pickle > None_exp_rate0.05349_scale0.0_notree_1_gammaFitAll_any_continuousF.notree.list
ls *notree_0.7_gammaFitAll*pickle > None_exp_rate0.05349_scale0.0_notree_0.7_gammaFitAll_any_continuousF.notree.list
ls *tree_1_exp_*_continuousF.run*pickle > None_exp_rate0.05349_scale0.0_tree_1_exp_any_continuousF.list
ls *tree_2.3_exp_*_continuousF.shortleaves*pickle > exp_node_order_rate_exp_rate0.1_scale0.0_tree_2.3_exp_any_continuousF.shortleaves.list
ls *longleaves*pickle > exp_node_order_rate_rev_exp_rate0.05_scale1.0tree_0.05_exp_any_continuousF.longleaves.list

#move to bucket
mv *list /bucket/KondrashovU/seq_space/simulate_evol/data/uniform_f_w_dist/
mv *csv /bucket/KondrashovU/seq_space/simulate_evol/data/uniform_f_w_dist/csv/ &
mv *fasta /bucket/KondrashovU/seq_space/simulate_evol/data/uniform_f_w_dist/fasta/ &
mv *pickle /bucket/KondrashovU/seq_space/simulate_evol/data/uniform_f_w_dist/trees/ &
mv *log /bucket/KondrashovU/seq_space/simulate_evol/data/uniform_f_w_dist/logs/ &
mv *pdistm /bucket/KondrashovU/seq_space/simulate_evol/data/uniform_f_w_dist/pdistm/ &

#calculate average distance from unique matrices
python ../../../scripts/get_ld_mean_median_var.py pdistm_unique_renamed/ simu_w_anc_mean_median_var_dist_unique &> simu_w_anc_get_ld_mean_median_var_pdistm_unique.log &
mv allstats_mean_median_var_ld_simu_w_anc_mean_median_var_dist_unique.csv csv/

#calculate dim from unique matrices
sbatch  get_dim_unique_03_range_from_data_w_anc.sh
cd /flash/KondrashovU/ladaisa/simulate_evol/uniform_frac_w_dist
mv kcoefs_* /bucket/KondrashovU/seq_space/simulate_evol/data/uniform_f_w_dist/
awk 'FNR==1 && NR!=1 {next} {print}' kcoefs_unique/*.coefnr > simu_w_anc_unique_dim_winsize0.3_12444.coefnr
awk 'FNR==1 && NR!=1 {next} {print}' kcoefs_range_from_data_unique/*.coefnr > simu_w_anc_unique_dim_winsize0.5_dim_range_from_data_12444.coefnr

#calculate dim with regression range from the data
#(this also calculates dnds from tree, but it's slow, so it's parallelized more below in a separate run)
# win_size = 0.5
sbatch get_dim_range_from_data_and_dnds_from_tree.sh
# win_size = 0.3 & 0.4 
sbatch ../../get_dim_unique_range_from_data_0304_w_anc.sh

#move to bucket
cd /flash/KondrashovU/ladaisa/simulate_evol/uniform_frac_w_dist/
mv *winsize0.5_dim_range_from_data.csv /bucket/KondrashovU/seq_space/simulate_evol/data/uniform_f_w_dist/csv/

mv /flash/KondrashovU/ladaisa/simulate_evol/uniform_frac_w_dist/kcoefs_range_from_data_unique/*0.4* kcoefs_range_from_data_unique_0.4/
mv /flash/KondrashovU/ladaisa/simulate_evol/uniform_frac_w_dist/kcoefs_range_from_data_unique/*0.3* kcoefs_range_from_data_unique_0.3/

#concatenate
awk 'FNR==1 && NR!=1 {next} {print}' *0.5_dim_range_from_data.csv > dist_dim_winsize0.5_dim_range_from_data_5244.csv
cp dist_dim_winsize0.5_dim_range_from_data_5244.csv ../../allstats/
# 0.3 & 0.4
cd kcoefs_range_from_data_unique_0.3/
awk 'FNR==1 && NR!=1 {next} {print}' *.coefnr > ../simu_kcoefs_mult_ali_max_k_range_from_data_non_log_0.3_20.coefnr
cd kcoefs_range_from_data_unique_0.4/
awk 'FNR==1 && NR!=1 {next} {print}' *.coefnr > ../simu_kcoefs_mult_ali_max_k_range_from_data_non_log_0.4_20.coefnr
cp simu_kcoefs_mult_ali_max_k_range_from_data_non_log_0.[3,4]_20.coefnr ../../../allstats/tmp/

#get dnds from tree
#split the list of trees first (script is slow)
split -l 50 --numeric-suffixes final_all_fastas_4845_notreefinal.list lists/final_all_fastas_4845_notreefinal.list
#calculate dnds from tree
sbatch get_dnds_from_tree_epistasis_notree_scaling_4845.sh

#move to bucket
cd /flash/KondrashovU/ladaisa/simulate_evol/uniform_frac_w_dist/
mv dnds_from_tree_list_final_trees_4845_notreefinal.*.tsv /bucket/KondrashovU/seq_space/simulate_evol/data/uniform_f_w_dist/csv/

#concatenate tsv
awk 'FNR==1 && NR!=1 {next} {print}' dnds_from_tree_list_final_trees_4845_notreefinal.*.tsv > dnds_from_tree_list_final_trees_4845_notreefinal_all.tsv
cp dnds_from_tree_list_final_trees_4845_notreefinal_all.tsv ../../allstats/

#keep original split tsv in a separate folder
cd /bucket/KondrashovU/seq_space/simulate_evol/data/uniform_f_w_dist/csv/
mkdir dnds_from_tree_split
mv dnds_from_tree_list_final_trees_4845_notreefinal.list*.tsv  dnds_from_tree_split/


#get dnds with PAML
sbatch  dnds_simulated_final_notree_tree_w_epi_paml.sh

#parse dnds
cd /flash/KondrashovU/ladaisa/simulate_evol/uniform_frac_w_dist/paml/ynout_dir/
ls *.ng86.ynout > /bucket/KondrashovU/seq_space/simulate_evol/data/uniform_f_w_dist/paml_dnds_ng86_simu.list
ls *.ynout |grep -v 'ng86'> /bucket/KondrashovU/seq_space/simulate_evol/data/uniform_f_w_dist/paml_dnds_yn00_simu.list

#create folders for output 
cd ../
mkdir ynout_dir_yn00_parsed_filtered
mkdir ynout_dir_ng86_parsed_filtered

#parse PAML outputs
sbatch  parse_dnds_yn00_ng86_simulated_frac_final.sh

#translate fastas into aa to get usage 
cd /bucket/KondrashovU/seq_space/simulate_evol/data/uniform_f_w_dist/fasta
for line in *fasta; do seqkit translate -T 1 $line > ../fasta_aa/$line; done 

#remove duplicates 
cd fasta_aa/
for line in *; do seqkit rmdup -s $line -o ../fasta_aa_unique/$line.no_duplicates; done  & 

#get usage from leaves
ml python/3.11.4
sbatch --partition compute -t 10:00:00 -c 1 --mem=20G --output=/flash/KondrashovU/ladaisa/logs/simulate_evol/get_orig_usage5244.log --job-name=simu_usage \
--wrap "python /bucket/KondrashovU/seq_space/scripts/get_usage_ali_seq_len_invar_sites_gap_frac.py  -f /bucket/KondrashovU/seq_space/simulate_evol/data/uniform_f_w_dist/final_all_fastas_5244.list -s '.fasta' -i /bucket/KondrashovU/seq_space/simulate_evol/data/uniform_f_w_dist/fasta_aa/ -o /flash/KondrashovU/ladaisa/simulate_evol/uniform_frac_w_dist/"

#unique (no duplicates)
sed 's/.fasta/.fasta.no_duplicates/g' final_all_fastas_5244.list  > final_all_fastas_5244_unique.list
sbatch --partition compute -t 10:00:00 -c 1 --mem=20G --output=/flash/KondrashovU/ladaisa/logs/simulate_evol/get_unique_usage5244.log --job-name=simu_usage_nd \
--wrap "python /bucket/KondrashovU/seq_space/scripts/get_usage_ali_seq_len_invar_sites_gap_frac.py  -f /bucket/KondrashovU/seq_space/simulate_evol/data/uniform_f_w_dist/final_all_fastas_5244_unique.list -s '.fasta.no_duplicates' -i /bucket/KondrashovU/seq_space/simulate_evol/data/uniform_f_w_dist/fasta_aa_unique/ -o /flash/KondrashovU/ladaisa/simulate_evol/uniform_frac_w_dist/"

#get usage from ancestors
ml python/3.11.4
sbatch --partition compute -t 10:00:00 -c 1 --mem=20G --output=/flash/KondrashovU/ladaisa/logs/simulate_evol/get_anc_usage12444.log --job-name=anc_usage \
--wrap "python /bucket/KondrashovU/seq_space/simulate_evol/get_usage_from_ancestors.py  -t /bucket/KondrashovU/seq_space/simulate_evol/data/uniform_f_w_dist/list_final_trees_5244.list -if /bucket/KondrashovU/seq_space/simulate_evol/data/uniform_f_w_dist/trees/ -of /flash/KondrashovU/ladaisa/simulate_evol/uniform_frac_w_dist/"



#get counts of the number of unique nonsyn substitutions (number of matrix samplings)
#(and number of syn subs and the total tree length)
cd /bucket/KondrashovU/seq_space/simulate_evol/
ml python/3.11.4
sbatch --partition compute -t 10:00:00 -c 1 --mem=20G --output=/flash/KondrashovU/ladaisa/logs/simulate_evol/trees5244_count_independent_nonsyn_syn_subs_treelen.log --job-name=simu_dnsum \
--wrap "python /bucket/KondrashovU/seq_space/simulate_evol/count_independent_nonsyn_syn_subs_treelen_from_trees.py -t /bucket/KondrashovU/seq_space/simulate_evol/data/uniform_f_w_dist/list_final_trees_5244.list -if /bucket/KondrashovU/seq_space/simulate_evol/data/uniform_f_w_dist/trees/ -of /flash/KondrashovU/ladaisa/simulate_evol/"


#get max pairwise distance
cd /bucket/KondrashovU/seq_space/simulate_evol/data/uniform_f_w_dist/pdistm_unique_renamed/
for matr in *pdistm; do awk -F',' '{for(i=1;i<=NF;i++)if($i>max)max=$i} END{sub(/\.pdistm$/, "", FILENAME); print FILENAME "," max}' ${matr} >> ../simu_max_dist_pdistm_matrices_unique_levenshtein.csv; echo $matr; done &> ../get_max_dist.log &


#---------------- Variable number of sequences simulations


#test dimensionality of the different number of randomly generated sequences from the matrix
#new matrix for each family 

ml python/3.11.4
sbatch --partition compute -t 10:00:00 -c 1 --array=100,1000 --mem=20G --output=/flash/KondrashovU/ladaisa/logs/simulate_evol/random_seq_fom_matrix_nseq10e2e3.log --job-name=nseqrand10e2e3 \
--wrap "python /bucket/KondrashovU/seq_space/simulate_evol/simulate_random_seqs_from_matrix_dim_for_diff_nseq.py -n \${SLURM_ARRAY_TASK_ID} -r 10 -w 0.3 -m /bucket/KondrashovU/seq_space/simulate_evol/params/fracs_and_matr_allowed_for_nseq_10linspace.pickle -of /flash/KondrashovU/ladaisa/simulate_evol/rand_seq_from_matr_nseqs/ &> /flash/KondrashovU/ladaisa/simulate_evol/random_seq_fom_matrix_nseq\${SLURM_ARRAY_TASK_ID}.log"

sbatch --partition compute -t 4-00:00:00 -c 1 --mem=100G --output=/flash/KondrashovU/ladaisa/logs/simulate_evol/random_seq_fom_matrix_nseq10e4.log --job-name=nseqrand10e4 \
--wrap "python /bucket/KondrashovU/seq_space/simulate_evol/simulate_random_seqs_from_matrix_dim_for_diff_nseq.py -n 10000 -r 10 -w 0.3 -m /bucket/KondrashovU/seq_space/simulate_evol/params/fracs_and_matr_allowed_for_nseq_10linspace.pickle -of /flash/KondrashovU/ladaisa/simulate_evol/rand_seq_from_matr_nseqs/ &> /flash/KondrashovU/ladaisa/simulate_evol/random_seq_fom_matrix_nseq_10000.log"

#concatenate 100, 1000 and 10000 result tables - read every wile ending with "csv", cat each line except for the first line of each file (excluding the first file)
awk 'FNR==1 && NR!=1 {next} {print}' *.csv > simu_random_sampling_nseq100-10000_repl10.csv
#and copy to allstats
cp  simu_random_sampling_nseq100-10000_repl10.csv  ../allstats/

sbatch --partition compute -t 4-00:00:00 -c 1 --mem=100G --output=/flash/KondrashovU/ladaisa/logs/simulate_evol/random_seq_fom_matrix_nseq50e4.log --job-name=nseqrand10e5 \
--wrap "python /bucket/KondrashovU/seq_space/simulate_evol/simulate_random_seqs_from_matrix_dim_for_diff_nseq.py -n 50000 -r 10 -w 0.3 -m /bucket/KondrashovU/seq_space/simulate_evol/params/fracs_and_matr_allowed_for_nseq_10linspace.pickle -of /flash/KondrashovU/ladaisa/simulate_evol/rand_seq_from_matr_nseqs/ &> /flash/KondrashovU/ladaisa/simulate_evol/random_seq_fom_matrix_nseq_50000.log"

sbatch --partition compute -t 4-00:00:00 -c 1 --mem=100G --output=/flash/KondrashovU/ladaisa/logs/simulate_evol/random_seq_fom_matrix_nseq10e5.log --job-name=nseqrand10e5 \
--wrap "python /bucket/KondrashovU/seq_space/simulate_evol/simulate_random_seqs_from_matrix_dim_for_diff_nseq.py -n 100000 -r 10 -w 0.3 -m /bucket/KondrashovU/seq_space/simulate_evol/params/fracs_and_matr_allowed_for_nseq_10linspace.pickle -of /flash/KondrashovU/ladaisa/simulate_evol/rand_seq_from_matr_nseqs/ &> /flash/KondrashovU/ladaisa/simulate_evol/random_seq_fom_matrix_nseq_100000.log"

#same matrix for each family
ml python/3.11.4
sbatch --partition compute -t 10:00:00 -c 1 --array=100,1000 --mem=20G --output=/flash/KondrashovU/ladaisa/logs/simulate_evol/random_seq_fom_matrix_nseq10e2e3.log --job-name=nseqrand10e2e3 \
--wrap "python /bucket/KondrashovU/seq_space/simulate_evol/simulate_random_seqs_from_matrix_dim_for_diff_nseq.py -n \${SLURM_ARRAY_TASK_ID} -r 10 -w 0.3 -nm False -m /bucket/KondrashovU/seq_space/simulate_evol/params/fracs_and_matr_allowed_for_nseq_10linspace.pickle -of /flash/KondrashovU/ladaisa/simulate_evol/rand_seq_from_matr_nseqs/ &> /flash/KondrashovU/ladaisa/simulate_evol/random_seq_fom_matrix_nseq\${SLURM_ARRAY_TASK_ID}.log"

sbatch --partition compute -t 4-00:00:00 -c 1 --mem=100G --output=/flash/KondrashovU/ladaisa/logs/simulate_evol/random_seq_fom_matrix_nseq10e4.log --job-name=nseqrand10e4 \
--wrap "python /bucket/KondrashovU/seq_space/simulate_evol/simulate_random_seqs_from_matrix_dim_for_diff_nseq.py -n 10000 -r 10 -w 0.3 -nm False -m /bucket/KondrashovU/seq_space/simulate_evol/params/fracs_and_matr_allowed_for_nseq_10linspace.pickle -of /flash/KondrashovU/ladaisa/simulate_evol/rand_seq_from_matr_nseqs/ &> /flash/KondrashovU/ladaisa/simulate_evol/random_seq_fom_matrix_nseq_10000.log"


#test dimensionality of the different number of sequences on notree
#same matrix for each family
ml python/3.11.4
sbatch --partition compute -t 4-00:00:00 -c 1 --array=100,1000 --mem=20G --output=/flash/KondrashovU/ladaisa/logs/simulate_evol/notree_seq_nseq10e2e3.log --job-name=nseqnotree10e2e3 \
--wrap "python /bucket/KondrashovU/seq_space/simulate_evol/simulate_tree_seqs_dim_for_diff_nseq.py -n \${SLURM_ARRAY_TASK_ID} -r 10 -w 0.3 -nm False -s True -m /bucket/KondrashovU/seq_space/simulate_evol/params/fracs_and_matr_allowed_for_nseq_10linspace.pickle -of /flash/KondrashovU/ladaisa/simulate_evol/rand_seq_from_matr_nseqs/ &> /flash/KondrashovU/ladaisa/simulate_evol/notree_seq_nseq\${SLURM_ARRAY_TASK_ID}.log"

sbatch --partition compute -t 4-00:00:00 -c 1 --mem=100G --output=/flash/KondrashovU/ladaisa/logs/simulate_evol/notree_seq_nseq10e4.log --job-name=nseqnotree10e4 \
--wrap "python /bucket/KondrashovU/seq_space/simulate_evol/simulate_tree_seqs_dim_for_diff_nseq.py -n 10000 -r 10 -w 0.3 -nm False -s True -m /bucket/KondrashovU/seq_space/simulate_evol/params/fracs_and_matr_allowed_for_nseq_10linspace.pickle -of /flash/KondrashovU/ladaisa/simulate_evol/rand_seq_from_matr_nseqs/ &> /flash/KondrashovU/ladaisa/simulate_evol/notree_seq_nseq_10000.log"



#test dimensionality of the different number of sequences on the tree
#same matrix for each family
ml python/3.11.4
sbatch --partition compute -t 4-00:00:00 -c 1 --array=100,1000 --mem=20G --output=/flash/KondrashovU/ladaisa/logs/simulate_evol/tree_seq_nseq10e2e3.log --job-name=nseqtree10e2e3 \
--wrap "python /bucket/KondrashovU/seq_space/simulate_evol/simulate_tree_seqs_dim_for_diff_nseq.py -n \${SLURM_ARRAY_TASK_ID} -r 10 -w 0.3 -nm False -m /bucket/KondrashovU/seq_space/simulate_evol/params/fracs_and_matr_allowed_for_nseq_10linspace.pickle -of /flash/KondrashovU/ladaisa/simulate_evol/rand_seq_from_matr_nseqs/ &> /flash/KondrashovU/ladaisa/simulate_evol/tree_seq_nseq\${SLURM_ARRAY_TASK_ID}.log"

sbatch --partition compute -t 4-00:00:00 -c 1 --mem=100G --output=/flash/KondrashovU/ladaisa/logs/simulate_evol/tree_seq_nseq10e4.log --job-name=nseqtree10e4 \
--wrap "python /bucket/KondrashovU/seq_space/simulate_evol/simulate_tree_seqs_dim_for_diff_nseq.py -n 10000 -r 10 -w 0.3 -nm False -m /bucket/KondrashovU/seq_space/simulate_evol/params/fracs_and_matr_allowed_for_nseq_10linspace.pickle -of /flash/KondrashovU/ladaisa/simulate_evol/rand_seq_from_matr_nseqs/ &> /flash/KondrashovU/ladaisa/simulate_evol/tree_seq_nseq_10000.log"

mv rand_seq_from_matr_nseqs/* /bucket/KondrashovU/seq_space/simulate_evol/data/rand_seq_from_same_matr_nseqs/
cp rand_seq_from_same_matr_nseqs/*csv allstats/nseqs/


#subsample matrices from tree with 10000 sequences
ls *tree_f*10000*pdistm > simu_tree_nseqs10000_pdistm.list
ml python/3.11.4
sbatch --partition compute -t 4-00:00:00 -c 1 --mem=50G --output=/flash/KondrashovU/ladaisa/logs/simulate_evol/subsample10e4.log --job-name=subsample10e4 \
--wrap "python /bucket/KondrashovU/seq_space/simulate_evol/subsample_tree_seqs_dim_for_diff_nseq.py -w 0.3 -f /bucket/KondrashovU/seq_space/simulate_evol/data/rand_seq_from_same_matr_nseqs/simu_tree_nseqs10000_pdistm.list -if /bucket/KondrashovU/seq_space/simulate_evol/data/rand_seq_from_same_matr_nseqs/ -of /flash/KondrashovU/ladaisa/simulate_evol/rand_seq_from_matr_nseqs/ &> /flash/KondrashovU/ladaisa/simulate_evol/rand_seq_from_matr_nseqs/tree_subsample_nseq10000.log"

#subsample matrices from starlike tree with 10000 sequences
ls *tree_starlike_f*10000*pdistm > simu_tree_starlike_nseqs10000_pdistm.list
ml python/3.11.4
sbatch --partition compute -t 4-00:00:00 -c 1 --mem=50G --output=/flash/KondrashovU/ladaisa/logs/simulate_evol/subsample_starlike10e4.log --job-name=subsample_star10e4 \
--wrap "python /bucket/KondrashovU/seq_space/simulate_evol/subsample_tree_seqs_dim_for_diff_nseq.py -w 0.3 -f /bucket/KondrashovU/seq_space/simulate_evol/data/rand_seq_from_same_matr_nseqs/simu_tree_starlike_nseqs10000_pdistm.list -if /bucket/KondrashovU/seq_space/simulate_evol/data/rand_seq_from_same_matr_nseqs/ -of /flash/KondrashovU/ladaisa/simulate_evol/rand_seq_from_matr_nseqs/ &> /flash/KondrashovU/ladaisa/simulate_evol/rand_seq_from_matr_nseqs/notree_subsample_nseq10000.log"


#subsample and recalculate dim with reg range from data and winsize=0.5 for all Nseq (100,1000,10000)- for random seq, tree and notree
#random seq
ls simu_random_sampling*10000*pdistm > simu_random_sampling_nseqs10000_pdistm.list
ml python/3.11.4
sbatch --partition compute -t 4-00:00:00 -c 1 --mem=50G --output=/flash/KondrashovU/ladaisa/logs/simulate_evol/randseq_nseq_dim_ur.log --job-name=randseq_nseq_dim_ur \
--wrap "python /bucket/KondrashovU/seq_space/simulate_evol/subsample_tree_seqs_dim_for_diff_nseq.py -n 1 -w 0.5 -r True -s '[100,1000,10000]' -f /bucket/KondrashovU/seq_space/simulate_evol/data/rand_seq_from_same_matr_nseqs/simu_random_sampling_nseqs10000_pdistm.list -if /bucket/KondrashovU/seq_space/simulate_evol/data/rand_seq_from_same_matr_nseqs/ -of /flash/KondrashovU/ladaisa/simulate_evol/rand_seq_from_matr_nseqs/ &> /flash/KondrashovU/ladaisa/simulate_evol/rand_seq_from_matr_nseqs/randseq_nseq_dim_ur.log"

#starlike tree
ml python/3.11.4
sbatch --partition compute -t 4-00:00:00 -c 1 --mem=50G --output=/flash/KondrashovU/ladaisa/logs/simulate_evol/starlike_nseq_dim_ur.log --job-name=starlike_nseq_dim_ur \
--wrap "python /bucket/KondrashovU/seq_space/simulate_evol/subsample_tree_seqs_dim_for_diff_nseq.py -n 1 -w 0.5 -r True -s '[100,1000,10000]' -f /bucket/KondrashovU/seq_space/simulate_evol/data/rand_seq_from_same_matr_nseqs/simu_tree_starlike_nseqs10000_pdistm.list -if /bucket/KondrashovU/seq_space/simulate_evol/data/rand_seq_from_same_matr_nseqs/ -of /flash/KondrashovU/ladaisa/simulate_evol/rand_seq_from_matr_nseqs/ &> /flash/KondrashovU/ladaisa/simulate_evol/rand_seq_from_matr_nseqs/starlike_nseq_dim_ur.log"

#tree
ml python/3.11.4
sbatch --partition compute -t 4-00:00:00 -c 1 --mem=50G --output=/flash/KondrashovU/ladaisa/logs/simulate_evol/tree_nseq_dim_ur.log --job-name=tree_nseq_dim_ur \
--wrap "python /bucket/KondrashovU/seq_space/simulate_evol/subsample_tree_seqs_dim_for_diff_nseq.py -n 1 -w 0.5 -r True -s '[100,1000,10000]' -f /bucket/KondrashovU/seq_space/simulate_evol/data/rand_seq_from_same_matr_nseqs/simu_tree_nseqs10000_pdistm.list -if /bucket/KondrashovU/seq_space/simulate_evol/data/rand_seq_from_same_matr_nseqs/ -of /flash/KondrashovU/ladaisa/simulate_evol/rand_seq_from_matr_nseqs/ &> /flash/KondrashovU/ladaisa/simulate_evol/rand_seq_from_matr_nseqs/tree_nseq_dim_ur.log"


#get max pairwise distance from matrices with variable number of sequences
cd /bucket/KondrashovU/seq_space/simulate_evol/data/rand_seq_from_same_matr_nseqs
for matr in *nseq10000*pdistm; do awk -F',' '{for(i=1;i<=NF;i++)if($i>max)max=$i} END{sub(/\.pdistm$/, "", FILENAME); print FILENAME "," max}' ${matr} >> ../simu_max_dist_pdistm_matrices_nseqs10000.csv; echo $matr; done &> ../get_max_dist.log &



