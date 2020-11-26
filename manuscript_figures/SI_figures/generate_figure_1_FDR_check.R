# set the paths 
source("~/Covid-19/EnsembleMHC-Covid19/manuscript_figures/set_paths.R")
library(ggpubr)
library(ggplot2)
library(parallel)
library(data.table)
library(stringr)
library(wesanderson)
library(latex2exp)
library(reshape2)
library(patchwork)
library(dplyr)

min_norm <- function(x) {
  (max(x) - x) / (max(x) - min(x))
}
std <- function(x) {
  (x - mean(x)) / sd(x)
}

# calculate PPV
PPV <- function(x) {
  table(x)["target"] / sum(table(x))
}

# calculate F1
F1 <- function(x, num_of_pos) {
  PPV <- table(x)["target"] / sum(table(x))
  recall <- table(x)["target"] / num_of_pos
  2 * ((PPV * recall) / (PPV + recall))
}


# set algorithms vector
algos <- c("mhcflurry_affinity_percentile", "mhcflurry_presentation_score", "MixMHCpred", "netMHC_affinity", "netMHCpan_EL_affinity", "netstab_affinity", "pickpocket_affinity")

# load the P_sum matrix. This is the stored algorithm and allele specific score and FDR
load(paste0(data_path, "P_sum_median_1000_boot.R"))

# read in the coronavirus predictions
files <- list.files(path = tumor_data, pattern = "csv", full.names = T)

# calcualte percentile score mhcflurry presentation
xx <- fread(list.files(tumor_data, pattern = "pred.out", full.names = T))

# clean up allele names
xx$HLA <- str_remove(xx$HLA, "HLA-")

# calculate percentile score for MHCflurry 
mhcflurry_thres <- do.call(rbind, lapply(unique(xx$HLA), function(q) {
  tmp <- xx %>% dplyr::slice(which(HLA == q))
  c(q, quantile(tmp$mhcflurry_presentation_score, .98), quantile(tmp$mhcflurry_presentation_score, .995))
})) %>% data.frame()

# name columns
colnames(mhcflurry_thres) <- c("HLA", "2", "0.5")

# store the thresholds
mhcflurry_thres$`2` <- as.numeric(as.character(mhcflurry_thres$`2`))
mhcflurry_thres$`0.5` <- as.numeric(as.character(mhcflurry_thres$`0.5`))
######

# perfrom calculations for all samples
all_out <- do.call(rbind, lapply(files, function(f) {
  
  print(f)
  # read in the cell line data
  cell_line_data <- fread(f)
  cell_line_data$HLA <- str_remove(cell_line_data$HLA, "HLA-")

  true_counts <- cell_line_data %>%
    slice(-which(duplicated(peptide))) %>%
    pull(ident) %>%
    table()
  pos <- true_counts["target"]
  # calculate peptide probabilities
  # start by iterating through every unique HLA
  pep_probs <- lapply(unique(cell_line_data$HLA), function(w) {
    # create a tmp variable with that consists of the corona virus predictions for one allele
    tmp <- cell_line_data %>% dplyr::slice(which(cell_line_data$HLA == w)) %>% data.frame()
    # normalize the scores for the presentation score and pickpocket
    # both of these scores are not percentiles and the thresholds were bsaed on the normalized scores
    tmp$mhcflurry_presentation_score <- min_norm(tmp$mhcflurry_presentation_score)
    tmp$pickpocket_affinity <- min_norm(tmp$pickpocket_affinity)
    # coverent the gene name into a factor
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
    prob_combo <- do.call(rbind, prob_list) %>%
      data.frame() %>%
      group_by(peptide) %>%
      summarise(prob = prod(prob)) %>%
      merge(tmp[, c("peptide", "HLA", "ident")])
  })

  # bind all the pep_prob list generated in the previous step
  all_counts <- do.call(rbind, pep_probs)

  # make vector of all unique alleles
  sel_alleles <- unique(all_counts$HLA)


 
  F1_ensemble <- do.call(rbind, lapply(c(seq(.01,1,.01)), function(i) {
    prob_threshold <- i
    df_all_proteins <- all_counts %>%
      slice(which(prob <= prob_threshold)) %>%
      slice(-which(duplicated(peptide)))
    data.frame(
      algo = "EnsembleMHC",
      F1 = as.numeric(F1(df_all_proteins$ident, pos)),
      PPV = as.numeric(PPV(df_all_proteins$ident)),
      recall = as.numeric(table(df_all_proteins$ident)["target"] / pos),
      threshold = prob_threshold
    )
  }))


  df <- F1_ensemble
  df$cell_line <- str_extract(f, "(?<=\\/)(CL|GB|ME|OV).+(?=_combo)")
  df
}))


p1 <- all_out %>%
  ggplot(aes(y = F1, x = threshold, group = threshold)) +
  geom_boxplot() +
  scale_x_continuous(breaks = c(0, 0.05, .25, .5, .75, 1)) +
  theme(axis.text.x = element_text(colour = c("black", "red", "black", "black", "black", "black"))) +
  ggtitle("F1")

p2 <- all_out %>%
  ggplot(aes(y = PPV, x = threshold, group = threshold)) +
  geom_boxplot() +
  scale_x_continuous(breaks = c(0, 0.05, .25, .5, .75, 1)) +
  theme(axis.text.x = element_text(colour = c("black", "red", "black", "black", "black", "black"))) +
  ggtitle("PPV")


p3 <- all_out %>%
  ggplot(aes(y = recall, x = threshold, group = threshold)) +
  geom_boxplot() +
  scale_x_continuous(breaks = c(0, 0.05, .25, .5, .75, 1)) +
  theme(axis.text.x = element_text(colour = c("black", "red", "black", "black", "black", "black"))) +
  ggtitle("recall")


library(patchwork)

p1+p2+p3 + plot_layout(nrow = 3)

#ggsave(garg, filename = paste0(Ensemble_PATH, "/plots/main_figures/Figure_1.pdf"))



