# read in paths
source("~/Covid-19/EnsembleMHC-Covid19/manuscript_figures/set_paths.R")
library(dplyr)
library(data.table)
library(stringr)
min_norm <- function(x) {
  (max(x) - x) / (max(x) - min(x))
}
library(UpSetR)
library(parallel)
load(paste0(data_path, "P_sum_median_1000_boot.R"))
algos <- c("mhcflurry_affinity_percentile", "mhcflurry_presentation_score", "MixMHCpred", "netMHC_affinity", "netMHCpan_EL_affinity", "netstab_affinity", "pickpocket_affinity")


corona_data <- fread(paste0(data_path, "predicted_corona_peptides.csv"), stringsAsFactors = F)

# calculate peptide prob and shape data  into correct format
pep_probs <- lapply(unique(corona_data$HLA), function(w) {
  # create a tmp variable with that consists of the corona virus predictions for one allele
  tmp <- corona_data %>% dplyr::slice(which(corona_data$HLA == w)) %>%data.frame()
  # normalize the scores for the presentation score and pickpocket
  # both of these scores are not percentiles and the thresholds were bsaed on the normalized scores
  tmp$mhcflurry_presentation_score <- min_norm(tmp$mhcflurry_presentation_score)
  tmp$pickpocket_affinity <- min_norm(tmp$pickpocket_affinity)
  # coverent the gene name into a factor
  tmp$gene <- factor(tmp$gene)
  # print current HLA being processed
  print(w)
  # create the prob list which is the selection of peptides that fall within the score filter for one algorithm
  # assign the algorithm FDR to each peptide based on the benchmarking calculations
  prob_list <- lapply(colnames(tmp)[colnames(tmp) %in% unique(P_sum$algo)], function(q) {
    # Originally, the algorithms were assigned PPVs. These PPVs are converted to FDR through the relation PPV = 1 - FDR
    neg <- 1 - P_sum$PPV[which(P_sum$HLA == str_remove(w, pattern = "\\*") & P_sum$algo == q)]
    # the score threshold required for that algorithm at the HLA is recovered from P_sum matrix
    thres <- P_sum$value[which(P_sum$HLA == str_remove(w, pattern = "\\*") & P_sum$algo == q)]
    # select all peptides that fall within the scoring threshold for that algorithm at that allele
    peptides <- tmp$peptide[which(tmp[, q] <= thres)]
    # return a data frame consisting of selected peptides, teh algorithm FDR, and algorithm name
    data.frame(peptide = peptides, prob = neg, algo = q)
  })
  # combine all of the prob_lists in one dataframe, calculate the products of the FDRs associatied with detecting algorithms
  # merge with information regarding the gene and HLA
  prob_combo <- do.call(rbind, prob_list) %>% data.frame()
  peps <- unique(prob_combo$peptide)
  FDR_pass<- prob_combo%>%group_by(peptide) %>% summarise(prob=prod(prob)) %>%slice(which(prob<=.05)) %>%mutate(algo="EnsembleMHC") %>% select(peptide,algo)
  
  peptide_mat<-rbind(prob_combo%>%select(peptide,algo),FDR_pass)
  
  peptide_mat <- peptide_mat %>% slice(which(peptide%in%FDR_pass$peptide))
  
  tmp <- data.frame(peptide = peps, sapply(c("EnsembleMHC",algos), function(i) {
    as.numeric(peps%in%peptide_mat$peptide[which(peptide_mat$algo == i)])
  }))
  
})



# bind all the pep_prob list generated in the previous step
all_prob <- do.call(rbind, pep_probs)
all_prob$peptide <- as.character(all_prob$peptide)

# create upset plot
sel_up <- upset(all_prob,
  sets = algos, sets.bar.color = "#ACED00",
  order.by = "freq", empty.intersections = "on"
)

sel_up
