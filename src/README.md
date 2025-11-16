# Python binaries for any dataset

## Filter, align and find orthologs
- split_mmseq.py 
- filter_fasta_by_num_of_seq.py
- compare_nt_translations_to_aa.py
- compare_nt_translations_to_aa_eggnogs.py
- get_aa_corresponding_to_nt.py
- get_nt_for_gammaproteo_eggnogs.py - script to get nt sequences for EggNOGs from NCBI genome assemblies

## Alignment statistics
- get_usage_ali_seq_len_invar_sites_gap_frac.py
- seq_length_invar_sites_gap_content_from_fasta.py
- get_mean_median_seq_len_from_fasta.py
- get_ld_mean_median.py (older version)

## Pairwise distances calculation & handling
- hamming_distance_multiple_ali_universal.py
- hamming_distance_pairwise_ali_universal.py
- hamming_distance_wo_gaps.py
- check_pairwise_pdistm_for_repeated_dist.py
- filter_identical_seqs_from_pdistm.py
- get_pdist_distribution_plot.py
- get_ld_mean_median_var.py
- subsample_pdistm_matrix.py

## Correlation dimension
- dimension_steepest_curve_from_matrix.py

## dN/dS
- gap_stop_codons_for_paml.py
- dnds_pairwise_ndata_n_allinfile_universal.py
- parse_yn00_output.py

# simulate_evol/ - folder with binaries for simulations

- **Main script to generate protein families and get required statistics**
    - simulate_evol/simulate_protein_family_diff_tree_w_epistasis_evol_stats_and_ali.py

- **Get statistics from simulated trees**
    - simulate_evol/get_dnds_from_simulated_tree.py
    - simulate_evol/get_dim_dist_from_simulated_tree.py
    - simulate_evol/get_usage_from_ancestors.py
    - simulate_evol/count_independent_nonsyn_syn_subs_treelen_from_trees.py

- **Simulate different sample size**
    - simulate_evol/simulate_tree_seqs_dim_for_diff_nseq.py
    - simulate_evol/simulate_random_seqs_from_matrix_dim_for_diff_nseq.py
    - simulate_evol/subsample_tree_seqs_dim_for_diff_nseq.py

# connected_components_cluster/ - folder to quantify sequence clustering in sequence space as a function of epistasis (for simulated seqeunces)

- connected_components_cluster/integral_b0_simulate_evol.txt
- connected_components_cluster/norm_integral_b0_simulate_evol.txt
- connected_components_cluster/barcode_to_curve_b0.cpp
- connected_components_cluster/barcode_to_curve_b0
- connected_components_cluster/create_scripts_b0.cpp
- connected_components_cluster/create_scripts_b0
- connected_components_cluster/drawing_curve_b0_simulate_evol.txt
- connected_components_cluster/example/Result_b0curve_simulate_evol/SimuOG8862.b0curve
- connected_components_cluster/example/Gnuplot_script_b0_simulate_evol/SimuOG8862_script.p
- connected_components_cluster/example/Drawing_curve_b0_simulate_evol/SimuOG8862_drawing.png
- connected_components_cluster/example/Result_b0curve_simulate_evol/SimuOG10371.b0curve
- connected_components_cluster/example/Gnuplot_script_b0_simulate_evol/SimuOG10371_script.p
- connected_components_cluster/example/Drawing_curve_b0_simulate_evol/SimuOG10371_drawing.png

# seq_space_lib/ - package with functions to calculate dimensionality and measure parameters for this project

# simulate_prot_evol_lib/ - package with functions for simulation part of the project 
