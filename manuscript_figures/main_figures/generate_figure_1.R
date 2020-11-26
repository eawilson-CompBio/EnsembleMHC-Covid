# read in paths
source("~/Covid-19/EnsembleMHC-Covid19/manuscript_figures/set_paths.R")

dirname(sys.frame(1)$ofile)
library(ggpubr)
library(ggplot2)
library(parallel)
library(data.table)
library(stringr)
library(wesanderson)
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


# set algorithm name vector
algos <- c("mhcflurry_affinity_percentile", "mhcflurry_presentation_score", "MixMHCpred", "netMHC_affinity", "netMHCpan_EL_affinity", "netstab_affinity", "pickpocket_affinity")

# load the parameterization summary matrix. This is the stored algorithm and allele specific score and FDR
load(paste0(data_path, "P_sum_median_1000_boot.R"))

# read in the coronavirus predictions
files <- list.files(path = tumor_data, pattern = "csv", full.names = T)

# calculate percentile score mhcflurry presentation
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

# perform calculations for all samples
all_out <- do.call(rbind, lapply(files, function(f) {
  print(f)
  # read in the cell line data
  cell_line_data <- fread(f)
  cell_line_data$HLA <- str_remove(cell_line_data$HLA, "HLA-")

  # get top number of MS-observed peptides 
  true_counts <- cell_line_data %>%
    slice(-which(duplicated(peptide))) %>%
    pull(ident) %>%
    table()
  pos <- true_counts["target"]
  
  # calculate peptide probabilities
  # start by iterating through every unique HLA
  pep_probs <- lapply(unique(cell_line_data$HLA), function(w) {
    # create a tmp variable consisting of the corona virus predictions for one allele
    tmp <- cell_line_data %>%
      dplyr::slice(which(cell_line_data$HLA == w)) %>%
      data.frame()
    # normalize the scores for the presentation score and pickpocket
    # both of these scores are not percentiles and the thresholds were based on the normalized scores
    tmp$mhcflurry_presentation_score <- min_norm(tmp$mhcflurry_presentation_score)
    tmp$pickpocket_affinity <- min_norm(tmp$pickpocket_affinity)

    # print current HLA being processed
    print(w)
    # create the prob list which is the selection of peptides that fall within the score filter for one algorithm
    # assign the algorithm FDR to each peptide based on the benchmarking calculations
    prob_list <- lapply(colnames(tmp)[colnames(tmp) %in% unique(P_sum$algo)], function(q) {
      # Originally, the algorithms were assigned PPVs. These PPVs are converted to FDR through the relation PPV = 1 - FDR
      neg <- 1 - P_sum$PPV[which(P_sum$HLA == str_remove(w, pattern = "\\*") & P_sum$algo == q)]
      # the score threshold required for that algorithm at that alelle is recovered from P_sum matrix
      thres <- P_sum$value[which(P_sum$HLA == str_remove(w, pattern = "\\*") & P_sum$algo == q)]
      # select all peptides that fall within the scoring threshold for that algorithm at that allele
      peptides <- tmp$peptide[which(tmp[, q] <= thres)]
      # return a data frame consisting of selected peptides, the algorithm FDR, and algorithm name
      data.frame(peptide = peptides, prob = neg, algo = q)
    })
    # combine all of the prob_lists in one dataframe, calculate the products of the FDRs associatied with the detecting algorithms
    # merge with information regarding the gene and HLA
    prob_combo <- do.call(rbind, prob_list) %>%
      data.frame() %>%
      group_by(peptide) %>%
      summarise(prob = prod(prob)) %>%
      merge(tmp[, c("peptide", "HLA", "ident")])
  })

  # bind all elements of the pep_prob list generated in the previous step
  all_counts <- do.call(rbind, pep_probs)

  # make vector of all unique alleles
  sel_alleles <- unique(all_counts$HLA)


  # calclate PPV, recall, and F1 for algorithms with percentile score
  algo_percentile <- c("mhcflurry_affinity_percentile", "MixMHCpred", "netMHC_affinity", "netMHCpan_EL_affinity", "netstab_affinity")
  F1_out <- do.call(rbind, lapply(c(.5, 2), function(i) {
    thres <- i
    score <- sapply(algo_percentile, function(w) {
      tmp <- cell_line_data %>% data.frame()
      tmp <- tmp[, c("peptide", "HLA", w, "ident")]
      tmp <- tmp[which(tmp[, 3] <= thres), ]
      if (length(unique(tmp$peptide)) != nrow(tmp)) {
        tmp <- tmp %>% slice(-which(duplicated(peptide)))
      }
      c(f = as.numeric(F1(tmp$ident, pos)), per = as.numeric(PPV(tmp$ident)), rec = as.numeric(table(tmp$ident)["target"] / pos))
    })
    data.frame(algo = colnames(score), F1 = as.numeric(score["f", ]), PPV = as.numeric(score["per", ]), recall = as.numeric(score["rec", ]), threshold = thres)
  }))
  #give threshold qualitative names
  F1_out$qual <- "high"
  F1_out$qual[which(F1_out$threshold == .5)] <- "low"

  # calculate metrics for pickpocket affinity
  # convert pickpoct to ic 50
  cell_line_data$pickpocket_affinity <- 50000^(1 - cell_line_data$pickpocket_affinity)
  algo_BA <- c("pickpocket_affinity")
  F1_out_BA <- do.call(rbind, lapply(c(50, 500), function(i) {
    thres <- i
    score <- sapply(algo_BA, function(w) {
      tmp <- cell_line_data %>% data.frame()
      tmp <- tmp[, c("peptide", "HLA", w, "ident")]
      tmp <- tmp[which(tmp[, 3] <= thres), ]
      if (length(unique(tmp$peptide)) != nrow(tmp)) {
        tmp <- tmp %>% slice(-which(duplicated(peptide)))
      }
      c(f = as.numeric(F1(tmp$ident, pos)), per = as.numeric(PPV(tmp$ident)), rec = as.numeric(table(tmp$ident)["target"] / pos))
    })
    data.frame(algo = colnames(score), F1 = as.numeric(score["f", ]), PPV = as.numeric(score["per", ]), recall = as.numeric(score["rec", ]), threshold = thres)
  }))
  F1_out_BA$qual <- "high"
  F1_out_BA$qual[which(F1_out_BA$threshold == 50)] <- "low"


  # calculate metrics for MHCflurry presentation
  F1_out_MHC <- do.call(rbind, lapply(c(2, .5), function(i) {
    thres <- i
    tmp <- cell_line_data %>% data.frame()
    tmp <- tmp[, c("peptide", "HLA", "mhcflurry_presentation_score", "ident")]
    thres <- mhcflurry_thres[, c("HLA", as.character(i))] %>% slice(which(HLA %in% unique(tmp$HLA)))
    score <- thres[, as.character(i)]
    names(score) <- thres$HLA

    passed_thres <- do.call(rbind, lapply(unique(tmp$HLA), function(i) {
      tmp %>% slice(which(HLA == i & mhcflurry_presentation_score >= score[i]))
    })) %>%
      arrange(desc(mhcflurry_presentation_score)) %>%
      slice(-which(duplicated(peptide)))

    if (length(unique(passed_thres$peptide)) != nrow(passed_thres)) {
      passed_thres <- passed_thres %>% slice(-which(duplicated(peptide)))
    }
    score <- c(f = as.numeric(F1(passed_thres$ident, pos)), per = as.numeric(PPV(passed_thres$ident)), rec = as.numeric(table(passed_thres$ident)["target"] / pos))
    data.frame(algo = "mhcflurry_presentation_score", F1 = as.numeric(score["f"]), PPV = as.numeric(score["per"]), recall = as.numeric(score["rec"]), threshold = i)
  }))
  F1_out_MHC$qual <- "high"
  F1_out_MHC$qual[which(F1_out_MHC$threshold == 0.5)] <- "low"

  F1_ensemble <- do.call(rbind, lapply(c(.5, .05), function(i) {
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
  F1_ensemble$qual <- "high"
  F1_ensemble$qual[which(F1_ensemble$threshold == .05)] <- "low"


  df <- rbind(F1_out, F1_out_BA, F1_out_MHC, F1_ensemble)
  df$cell_line <- str_extract(f, "(?<=\\/)(CL|GB|ME|OV).+(?=_combo)")
  df
}))

# rename factors for plot
all_out$algo <- as.character(all_out$algo)
all_out$qual[which(all_out$qual == "low")] <- "restrictive"
all_out$qual[which(all_out$qual == "high")] <- "permissive"

# format algorithm names
all_out$algo[which(all_out$algo == "pickpocket_affinity")] <- "PickPocket"
all_out$algo[which(all_out$algo == "netstab_affinity")] <- "netMHCstabpan"
all_out$algo[which(all_out$algo == "netMHCpan_EL_affinity")] <- "netMHCpan-4.0-EL"
all_out$algo[which(all_out$algo == "netMHC_affinity")] <- "netMHC-4.0"
all_out$algo[which(all_out$algo == "mhcflurry_affinity_percentile")] <- "MHCflurry-affinity"
all_out$algo[which(all_out$algo == "mhcflurry_presentation_score")] <- "MHCflurry-presentation"

# format tumor data names
all_out$cell_line[which(all_out$cell_line == "CLL_DFCI-5283_2018")] <- "CLL C"
all_out$cell_line[which(all_out$cell_line == "CLL_DFCI-5328_20180512")] <- "CLL B"
all_out$cell_line[which(all_out$cell_line == "GBM_H4198_BT187")] <- "GBM 9"
all_out$cell_line[which(all_out$cell_line == "GBM_H4512_BT145")] <- "GBM 11"
all_out$cell_line[which(all_out$cell_line == "GBM_Pat7")] <- "GBM 7"
all_out$cell_line[which(all_out$cell_line == "MEL_13240-002")] <- "MEL 1"
all_out$cell_line[which(all_out$cell_line == "MEL_13240-005")] <- "MEL 2"
all_out$cell_line[which(all_out$cell_line == "MEL_13240-006")] <- "MEL 3"
all_out$cell_line[which(all_out$cell_line == "MEL_13240-015")] <- "MEL 15"
all_out$cell_line[which(all_out$cell_line == "OV_CP-594_v1_20161007")] <- "OV 1"

#----------------------------------------
# generate plots
#----------------------------------------


all_out$combo_name <- paste(all_out$algo, all_out$qual, sep = "_")
lev <- all_out$combo_name[which(all_out$cell_line == "CLL C")]
all_out$combo_name <- factor(all_out$combo_name, levels = rev(lev))
cols <- wes_palette(name = "Zissou1", type = "continuous")
lab_col <- c(rep("#663333", 2), rep("#ff6600", 2), rep("#006600", 2), rep("#3399cc", 10))
all_out$cell_line <- factor(all_out$cell_line, unique(all_out$cell_line)[c(2, 1, 5, 3, 4, 6, 7, 8, 9, 10)])

colnames(all_out)[3] <- "precision"


p_b_bar <- all_out %>%
  dplyr::select(algo, precision, recall, threshold, qual) %>%
  group_by(algo, qual, threshold) %>%
  mutate(precision = mean(precision, na.rm = T), recall = mean(recall, na.rm = T)) %>%
  melt(id.vars = c("algo", "qual", "threshold")) %>%
  unique() %>%
  ggplot(aes(x = algo, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap(qual ~ ., ncol = 1) +
  labs(fill = "metric") +
  xlab("") +
  ylab("") +
  theme_pubclean() +
  theme(axis.text.x = element_text(angle = 90)) +
  guides(fill = guide_legend(title.position = "top", title.hjust = .5))


hp <- all_out %>% ggplot(aes(y = combo_name, x = cell_line, fill = F1)) +
  geom_tile() +
  theme_classic() +
  scale_fill_gradientn(colours = cols) +
  theme(legend.position = "top", axis.text.x = element_text(angle = 90), axis.text.y = element_text(colour = lab_col)) +
  xlab("") +
  ylab("") +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = .5))

bp <- all_out %>%
  slice(-which(is.na(F1))) %>%
  group_by(combo_name) %>%
  summarise(F1 = mean(F1)) %>%
  ggplot(aes(y = combo_name, x = F1)) +
  geom_bar(stat = "identity") +
  theme_pubclean() +
  theme(axis.text.y = element_blank()) +
  ylab("") +
  xlab("precision")

garg <- p_b_bar + {
  hp + bp + plot_layout(widths = c(.7, .3))
}

garg

