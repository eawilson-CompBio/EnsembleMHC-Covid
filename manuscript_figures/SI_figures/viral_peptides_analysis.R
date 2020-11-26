# read in paths
source("~/Covid-19/EnsembleMHC-Covid19/manuscript_figures/set_paths.R")
library(stringr)
library(dplyr)
library(parallel)
library(ggplot2)
library(data.table)
library(RColorBrewer)
library(ggthemr)

# functions ---------------------------------------------------------------

min_norm <- function(x) {
  (max(x) - x) / (max(x) - min(x))
}

eval_metrics <- function(ref_peps, target_peps) {
  hits <- sum(sapply(unique(ref_peps$MHC), function(i) {
    tmp_ref <- ref_peps$peptide[which(ref_peps$MHC == i)]
    tmp_target <- target_peps$peptide[which(target_peps$HLA == i)]
    sum(tmp_target %in% tmp_ref)
  }))
  
  recall <- hits / length(ref_peps$peptide)
  precision <- hits / length(target_peps$peptide)
  
  F1 <- 2 * ((precision * recall) / (precision + recall))
  data.frame(precision, recall, F1)
}


netCTL <- fread(paste0(viral_peptide_dir,"/all_out_netCTL_trimmed.txt"),header = F,sep = " ")
netMHCcons <- fread(paste0(viral_peptide_dir,"/netMHCcons_all_out.txt"),header = F,sep = " ")
colnames(netCTL) <-c("protein","HLA","peptide","affinity")
colnames(netMHCcons) <-c("HLA","peptide","affinity")
netCTL$HLA<-str_remove_all(netCTL$HLA,"HLA-|\\*")
netMHCcons$HLA<-str_remove_all(netMHCcons$HLA,"HLA-|\\*")

select_peptides<-read.csv(paste0(viral_peptide_dir,"/immunogenic_viral_peptides.csv"))

viral_pred_peptides <- do.call(rbind, mclapply(list.files(paste0(viral_peptide_dir), pattern = "_scored_peptides.csv", full.names = T), function(i) {
  data.table::fread(i)
}, mc.cores = 10)) %>%
  data.frame() %>%
  mutate(pickpocket_affinity = 50000^(1 - pickpocket_affinity))

viral_pred_peptides$HLA <- str_remove_all(viral_pred_peptides$HLA, "HLA-|\\*")
load(paste0(data_path, "P_sum_median_1000_boot.R"))


# only look at peptides that there is some HLA data

sel_viral_HLA <- unique(select_peptides$MHC)

pep_probs <- lapply(sel_viral_HLA, function(w) {
  # create a tmp variable with that consists of the corona virus predictions for one allele
  tmp <- viral_pred_peptides %>% dplyr::slice(which(viral_pred_peptides$HLA == w))
  # normalize the scores for the presentation score and pickpocket
  # both of these scores are not percentiles and the identifed thresholds were based on the similarly normalized scores
  tmp$mhcflurry_presentation_score <- min_norm(tmp$mhcflurry_presentation_score)
  tmp$pickpocket_affinity <- min_norm(tmp$pickpocket_affinity)
  # coverent the gene name into a factor
  tmp$gene <- factor(tmp$gene)
  # print current HLA being processed
  print(w)
  # create the prob list which is the selection of peptides that fall within the score filter for one algorithm
  # assign the algorithm FDR to each peptide based on the benchmarking calculations
  # loop through all 7 algorithms
  prob_list <- lapply(colnames(tmp)[colnames(tmp) %in% unique(P_sum$algo)], function(q) {
    # Originally, the algorithms were assigned PPVs. These PPVs are converted to FDR through the relation FDR = 1 - PPV
    # therefore, the following line finds the PPV for the allele specified by w and algorithm q
    neg <- 1 - P_sum$PPV[which(P_sum$HLA == str_remove(w, pattern = "\\*") & P_sum$algo == q)]
    # same as above but with the score threshold for teh w and q combo
    thres <- P_sum$value[which(P_sum$HLA == str_remove(w, pattern = "\\*") & P_sum$algo == q)]
    # select all peptides that fall within the scoring threshold for that algorithm at that allele
    peptides <- tmp$peptide[which(tmp[, q] <= thres)]
    # return a data frame consisting of selected peptides, the algorithm FDR, and algorithm name
    data.frame(peptide = peptides, prob = neg, algo = q)
  })
  # combine all of the prob_list elements in one dataframe, calculate the products of the FDRs associatied with detecting algorithms
  # merge with information regarding the gene and HLA
  prob_combo <- do.call(rbind, prob_list) %>%
    data.frame() %>%
    group_by(peptide) %>%
    summarise(prob = prod(prob)) %>%
    merge(tmp[, c("peptide", "gene", "HLA")])
})



# bind all the pep_prob list generated in the previous step
all_counts <- do.call(rbind, pep_probs)

prob_threshold <- .05

# Summarize by allele and filter for peptides that are less than or equal to 5%
# return a count of
pass_filter <- all_counts %>%
  slice(which(prob <= prob_threshold))




BA_algos <- colnames(viral_pred_peptides)[c(7)]
percentile_algos <- colnames(viral_pred_peptides)[c(2:3, 5:6, 12)]
mhcflurry_pres_algo <- "mhcflurry_presentation_score"

collected_low_high <- do.call(rbind, mclapply(seq(25, nrow(select_peptides), 10), function(num_peps) {
  low_to_high_algos <- do.call(rbind, lapply(c(BA_algos, percentile_algos), function(i) {
    target <- viral_pred_peptides[, c("peptide", "HLA", i)]
    target <- target[order(target[[i]]), ][1:num_peps, ]
    eval_metrics(ref_peps = select_peptides, target_peps = target) %>% data.frame(., algo = i, num_peps = num_peps)
  }))
}, mc.cores = 10))

collected_low_high_mhcflurry <- do.call(rbind, mclapply(seq(25, nrow(select_peptides), 10), function(num_peps) {
  target <- viral_pred_peptides[, c("peptide", "HLA", "mhcflurry_presentation_score")]
  target <- target[order(target$mhcflurry_presentation_score, decreasing = T), ][1:num_peps, ]
  mhcflurry_top <- eval_metrics(ref_peps = select_peptides, target_peps = target) %>% data.frame(., algo = "mhcflurry_presentation_score", num_peps = num_peps)
}, mc.cores = 10))

collected_ensemble <- do.call(rbind, mclapply(seq(25, nrow(select_peptides), 10), function(num_peps) {
  target <- all_counts %>%
    arrange(prob) %>%
    slice(1:num_peps)
  Ensemble_top <- eval_metrics(ref_peps = select_peptides, target_peps = target) %>% data.frame(., algo = "EnsembleMHC", num_peps = num_peps)
}, mc.cores = 10))



collected_netCTL <- do.call(rbind, mclapply(seq(25, nrow(select_peptides), 10), function(num_peps) {
   target <- netCTL[, c("peptide", "HLA","affinity")]
   target <- target[order(target$affinity), ][1:num_peps, ]
   eval_metrics(ref_peps = select_peptides, target_peps = target) %>% data.frame(., algo = "netCTLpan", num_peps = num_peps)

}, mc.cores = 10))

collected_netMHCcons <- do.call(rbind, mclapply(seq(25, nrow(select_peptides), 10), function(num_peps) {
  target <- netMHCcons[, c("peptide", "HLA","affinity")]
  target <- target[order(target$affinity), ][1:num_peps, ]
  eval_metrics(ref_peps = select_peptides, target_peps = target) %>% data.frame(., algo = "netMHCcons", num_peps = num_peps)
  
}, mc.cores = 10))


df_top_pep <- rbind(collected_low_high, 
                    collected_low_high_mhcflurry, 
                    collected_ensemble, 
                    collected_netCTL,
                    collected_netMHCcons)


ggthemr("fresh")
df_top_pep %>%
  reshape2::melt(id.vars = c("algo", "num_peps")) %>%
  slice(which(variable == "precision")) %>%
  ggplot(aes(x = num_peps, y = as.numeric(value), color = algo)) +
  geom_line() +
  scale_color_manual(values = brewer.pal(10, "Paired")) +
  theme(
    legend.position = "top",
    legend.title = element_blank()
  ) +
  guides(color = guide_legend(override.aes = list(size = 7))) +
  ylab("precision") +
  xlab("number of peptides (ranked by affinity)")


