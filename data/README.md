# Data

This folder contains summary data tables which can be used to replicate plots.

## ./ - main datasets of natural and simulated protein families

- vertebrates_summary.csv
- enterobacterales_summary.csv
- gammaproteobacteria_summary.csv
- cogs_summary.csv
- simulated_summary.csv

### List of columns in the datasets:

- **ID columns:**
    - OG_id
    - SimuID (simulated OGs only)

- **Alignment/OG statistics:**
    - Number of sequences
        - num_of_seq_aa_orig 
        - num_of_seq_aa_unique
    - Alignemnt and sequence length
        - ali_len
        - mean_seq_len_orig
        - median_seq_len_orig
    - Variable and invariant sites
        - num_invar_sites
        - num_invar_sites_less02gaps
        - num_invar_sites_less05gaps
        - num_invar_sites_no_gaps
        - num_var_sites
        - num_var_sites_less02gaps
        - num_var_sites_less05gaps

    - Percentage of gaps
        - num_sites_less02gaps
        - num_sites_less05gaps
        - median_frac_gaps_per_site
        - mean_frac_gaps_per_site


- **Amino acid usage:**
    - mean_alphabet_size
    - mean_alphabet_size_less02gaps
    - mean_alphabet_size_less05gaps
    - median_alphabet_size
    - median_alphabet_size_less02gaps
    - median_alphabet_size_less05gaps

- **Levenshtein Distance matrix statistics:**
    - max_ld
    - mean_ld
    - median_ld
    - var_ld

- **Correlation dimension:**
    - dim_ur
    - dim_CI_ur
    - dim_R^2_ur
    - dim_range_end_ur
    - dim_range_start_ur

- **Expected and observed effective topological dimension**
    - Observed:
        - n_seq_from_dim_log10
        - n_seq_from_dim_log20
        - n_seq_from_dim_log_alpha
        - log10_nseq_ci
        - log20_nseq_ci

    - Expected
        - n_seq_bin_expected_from_dn/ds_log20
        - n_seq_bin_expected_from_usage_log20
        - n_seq_expected_from_dn/ds_log20
        - n_seq_expected_from_usage_log20

    - Difference between expected and observed dimension:
        - n_seq_diff_dnds_log10
        - n_seq_diff_dnds_log20
        - n_seq_diff_usage_log10
        - n_seq_diff_usage_log20


- **Observed dN/dS from PAML**
    - dnds_mean
    - dnds_median
    - dnds_var
    - dn_mean
    - dn_median
    - dn_var
    - ds_mean
    - ds_median
    - ds_var
    - filt_length
    - frac_NA
    - frac_SE_high
    - frac_dS<0.1
    - frac_dS>0.8
    - orig_length

- **Expected dN/dS calculated from amino acid usage**
    - dnds_from_usage
    - dnds_from_usage_cons
    - dnds_from_usage_sites80
    - dnds_from_usage_sites80_cons

- **Clustering measure (simulated only)**
    - integral_b0
    - norm_integral_b0

- **Annotation information**
    - Prokaryotic datasets - inherited from corresponding COG
        - COG
        - func_gr_all
        - func_gr
        - descr
    - Vertebrates
        - geneID

- **Simulation parameters:**
    - Inputs
        - nodes
        - branch_length_distr
        - fraction_allowed_subs
        - gamma
        - mode
        - replicate
        - run
        - suffix
        - tree_len_scale
    - Measured statistics from tree: N, S, dN/dS, N+S 
        - dnds_mean_tree
        - dnds_median_tree
        - dnds_var_tree
        - nonsyn_subs_count
        - syn_subs_count
        - total_treelength



## small_simulations/ folder - small datasets simulated to test the effect of different parameters (sequence length, amino acid usage, sample size etc.) on dimensionality estimates

### Dimensionality of simulated random fit sequences and sequence from starlike and binary tree for alpha:[0.07, 1] 

- small_simulations/simu_compare_random_trees_F0.1-1_rep10_nseq1000.csv
- small_simulations/simu_compare_random_trees_F0.07-0.6_rep10_nseq1000.csv

### Dimension as a function of sequence length and amino acid usage 

- small_simulations/random_seq_dim_length2-15_usage2-20_map.csv
- small_simulations/random_seq_dim_length50-500_usage2-20_map.csv
- small_simulations/random_seq_from_aa_dim_winsize03-05_length100-500_usage2-20_map.csv

### Effect of the sample size for dimension estimation

- small_simulations/simu_random_sampling_nseq100_1000_10000_subsampled_from10000_dim_range_from_data.csv
- small_simulations/simu_tree_nseq100_1000_10000_subsampled_from10000_dim_range_from_data.csv
- small_simulations/simu_tree_starlike_nseq100_1000_10000_subsampled_from10000_dim_range_from_data.csv
- small_simulations/simu_max_dist_pdistm_matrices_nseqs10000.csv 

### Sequences evolving on a starlike tree for increasingly longer amounts of time approach dimensionality of randomly sampled sequences

- small_simulations/rand_and_starlike_long_evol_x2-10_w_syn_counts.csv

## dim_diff_reg_range/ folder - testing alternative ways to estimate correlation dimension

### Different regression window size

- dim_diff_reg_range/eggnog_gamma_kcoefs_mult_ali_max_k_range_from_data_non_log_0.3_20.coefnr
- dim_diff_reg_range/eggnog_gamma_kcoefs_mult_ali_max_k_range_from_data_non_log_0.4_20.coefnr
- dim_diff_reg_range/entero_kcoefs_mult_ali_max_k_range_from_data_non_log_0.3_20.coefnr
- dim_diff_reg_range/entero_kcoefs_mult_ali_max_k_range_from_data_non_log_0.4_20.coefnr
- dim_diff_reg_range/kcoefs_mult_ali_max_k_range_from_data_unique_non_log_0.3_20.coefnr
- dim_diff_reg_range/kcoefs_mult_ali_max_k_range_from_data_unique_non_log_0.4_20.coefnr
- dim_diff_reg_range/simu_kcoefs_mult_ali_max_k_range_from_data_non_log_0.3_20.coefnr
- dim_diff_reg_range/simu_kcoefs_mult_ali_max_k_range_from_data_non_log_0.4_20.coefnr

### Distances from pairwise alignemnts

- dim_diff_reg_range/entero_pairali_kcoefs_max_k_range_from_data_non_log_0.5_20.coefnr
- dim_diff_reg_range/gamma_pairali_kcoefs_max_k_range_from_data_non_log_0.5_20.coefnr
- dim_diff_reg_range/vert_pairali_kcoefs_max_k_range_from_data_non_log_0.5_20.coefnr

- dim_diff_reg_range/random_0.05_og_ids_entero_for_pairali.list
- dim_diff_reg_range/random_0.05_og_ids_gamma_for_pairali.list
- dim_diff_reg_range/random_0.05_og_ids_vert_for_pairali.list

## examples/ folder - sample pairwise distance matrices

## distances/ folder - sample pairwise distance distributions from each orthogroup in each dataset

## metadata/ folder - additional information about datasets (genome NCBI IDs, etc.)

- metadata/Enterobacterales_genomes501.tsv

## Raw data

Raw data can be found at Zenodo DOI:10.5281/zenodo.17623392

Raw data contains:

**Vertebrates dataset - ncbi_vertebrates_data.zip**
- Amino acid alignemnts (created as a part of this project)
- Nucleotide alignemnts for dN/dS calculation (created as a part of this project)
- Raw amino acid sequences of sample 616 protein families used for pairwise alignments 
(unlike other datasets can't be obtained by removing the gaps from alignments due to use of 
MACSE for codon alignment)
- Species annotation file for sequence ids 

**Enterobacterales dataset - ncbi_enterobacterales_data.zip**
- Amino acid alignemnts (created as a part of this project)
- Nucleotide alignemnts for dN/dS calculation (created as a part of this project)
- Species annotation file for sequence ids 

**Gammaproteobacteria dataset - eggnog_gammaproteobacteria.zip**
- Nucleotide alignemnts for dN/dS calculation (created as a part of this project)

**COGs**
- Pairwise distance matrices (sequences where pairwise-aligned "on the run" while distances were 
calculated, see ../src/hamming_distance_pairwise_ali_universal.py)


Other data not provided (pairwise distances, dN/dS, alignment statistics) can be calculated from
provided/publicly available alignments (as is the case for Gammaproteobacteria EggNOG alignements) -
see ../workflows/pipeline.sh for more details.