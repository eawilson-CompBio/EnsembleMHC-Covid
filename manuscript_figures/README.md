# Manuscript Figure generation

This directory is primarily for recreating the manuscript figures. As a general rule of thumb, generate the main text figures before the SI figures. 

## setup
 
 In order for the scripts to run properly, paths to necessary datasets need to be set.
 
 * change the paths in the [set_path.R](set_path.R) script to link important datasets to each script.
 
 * modify the source call in the header of each script so that the set_path.R script is sourced. 
 
 
## datasets 

large datasets needed to replicate all figures can be downloaded from the following: https://drive.google.com/drive/folders/1FhA6wTpDfmnn2iO6eG09EYLGMfsYzhzK?usp=sharing



## Figure to script mapping

* Figure 1: generate_figure_1.R
* Figure 2: generate_figure_2.R
* Figure 3: generate_figure_3.R
* Figure 4: generate_figure_4.R
* SI A2: SI_figure_algotihm_benchmarking.R
* SI A3: generate_figure_1_FDR_check.R
* SI A4: generate_figure_1_MCC.R
* SI A5: generate_figure_1_with_netMHCcons.R & Viral_peptides_analysis.R
* SI A7: SI_upset_plots.R
* SI A8: SI_pep_length_dist.R 
* SI A9: SI_generate_logo_plots.R
* SI A11: SI_analyze_mutants.R
* SI A12: SI_even_dist.R
* SI A13: peptide_overlap_withother_papers.R
* SI A14: SI_structural_and_all_proteins_corr.R
* SI A15-16: SI_CI.R
* SI A17/25: SI_alt_corr_methods_pwr.R
* SI A18: test_downsampling_figure_3.R
* SI A19: SI_HLA_shuffle.R 
* SI A20: SI_individual_algorithms.R
* SI A21: SI_compare_bin_and_single_lm.R
* SI A22: SI_diff_allele_normalization.R
* SI A24: SI_generate_qq.R
