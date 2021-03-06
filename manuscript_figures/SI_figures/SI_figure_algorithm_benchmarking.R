source("~/Covid-19/EnsembleMHC-Covid19/manuscript_figures/set_paths.R")
library(ggthemes)
library(ggpubr)
library(ggplot2)
library(dplyr)
library(data.table)
library(stringr)
library(parallel)
library(wesanderson)
min_norm <- function(x) {
  (max(x) - x) / (max(x) - min(x))
}
std <- function(x) {
  (x - mean(x)) / sd(x)
}
library(latex2exp)
library(ggridges)
library(ggrepel)
library(reshape2)

# set the selected algorithms that were available for each prediction
algos <- c("mhcflurry_affinity_percentile", "mhcflurry_presentation_score", "MixMHCpred", "netMHC_affinity", "netMHCpan_EL_affinity", "netstab_affinity", "pickpocket_affinity")
# read in the benchmariking files
benchmark_files <- list.files(path = paste0(data_path, "allele_benchmark_datasets/"), pattern = "csv", full.names = T)
# load in vetor with all selected alleles
load(paste0(data_path, "selected_alleles.R"))
# load in the P_sum matrix with median value for PPV and score threshold
load(paste0(data_path, "P_sum_median_1000_boot.R"))

# get path to all of the benchmarking files

# create correlation plot
all <- mclapply(benchmark_files, function(file) {
  # read in benchmark with index file
  mat <- read.csv(file, stringsAsFactors = F)
  # make sure there are no duplicated rrrow
  mat <- mat %>% unique()
  # normailze score for mhcflurry
  mat$mhcflurry_presentation_score <- min_norm(mat$mhcflurry_presentation_score)
  # normalize score for pickpocket
  mat$pickpocket_affinity <- min_norm(mat$pickpocket_affinity)
  # return matrix
  mat
}, mc.cores = 10)

# find rf files with at least 20 column inidcating prediction for all alleles
all_sel <- do.call(rbind, all[which(sapply(all, ncol) == 20)])
# sort by peptide then HLA to ensure correct arrangement
all_sel <- all_sel %>%
  arrange(peptide) %>%
  arrange(HLA)
# calculate score correlation
cor_all <- cor(all_sel[, which(colnames(all_sel) %in% algos)])

# rename the algorithm to pretty names
P_sum$algo[which(P_sum$algo == "pickpocket_affinity")] <- "PickPocket"
P_sum$algo[which(P_sum$algo == "netstab_affinity")] <- "netMHCstabpan"
P_sum$algo[which(P_sum$algo == "netMHCpan_EL_affinity")] <- "netMHCpan-4.0-EL"
P_sum$algo[which(P_sum$algo == "netMHC_affinity")] <- "netMHC-4.0"
P_sum$algo[which(P_sum$algo == "mhcflurry_affinity_percentile")] <- "MHCflurry-affinity"
P_sum$algo[which(P_sum$algo == "mhcflurry_presentation_score")] <- "MHCflurry-presentation"

rownames(cor_all)[which(rownames(cor_all) == "pickpocket_affinity")] <- "PickPocket"
rownames(cor_all)[which(rownames(cor_all) == "netstab_affinity")] <- "netMHCstabpan"
rownames(cor_all)[which(rownames(cor_all) == "netMHCpan_EL_affinity")] <- "netMHCpan-4.0-EL"
rownames(cor_all)[which(rownames(cor_all) == "netMHC_affinity")] <- "netMHC-4.0"
rownames(cor_all)[which(rownames(cor_all) == "mhcflurry_affinity_percentile")] <- "MHCflurry-affinity"
rownames(cor_all)[which(rownames(cor_all) == "mhcflurry_presentation_score")] <- "MHCflurry-presentation"

colnames(cor_all)[which(colnames(cor_all) == "pickpocket_affinity")] <- "PickPocket"
colnames(cor_all)[which(colnames(cor_all) == "netstab_affinity")] <- "netMHCstabpan"
colnames(cor_all)[which(colnames(cor_all) == "netMHCpan_EL_affinity")] <- "netMHCpan-4.0-EL"
colnames(cor_all)[which(colnames(cor_all) == "netMHC_affinity")] <- "netMHC-4.0"
colnames(cor_all)[which(colnames(cor_all) == "mhcflurry_affinity_percentile")] <- "MHCflurry-affinity"
colnames(cor_all)[which(colnames(cor_all) == "mhcflurry_presentation_score")] <- "MHCflurry-presentation"


# set colors for correlation matrix
cols <- wes_palette(name = "Zissou1", type = "continuous")
# make correlation matrix figure
cor_plot <- melt(cor_all) %>% ggplot(aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradientn(colours = cols) +
  geom_text(aes(label = round(value, 2))) +
  theme_pubclean() +
  theme(axis.text.x = element_text(angle = 90)) +
  xlab("") +
  ylab("") +
  labs(fill = "correlation")

# calculate the per algorithm FDR
# note: scores are orginally reproted as the PPV, and are transfromed to FDR through the relation of FDR = 1 - PPV
algo_FDR <- P_sum %>% ggplot(aes(x = 1 - PPV, y = algo, fill = algo)) +
  geom_density_ridges() +
  ylab("algorithm") +
  theme(legend.position = "none") +
  xlab("False Detection Rate")

# calculate the per allele FDR
by_allele_FDR <- P_sum %>%
  slice(which(HLA %in% sel_alleles)) %>%
  ggplot(aes(x = factor(HLA, levels = rev(unique(HLA))), y = 1 - PPV)) +
  geom_boxplot(fill = "#00ACED") +
  coord_flip() +
  theme_classic() +
  ylab("False Detection Rate") +
  xlab("HLA")

# arrange plots
SI_fig <- ggarrange(
  by_allele_FDR,
  ggarrange(
    algo_FDR + theme_classic() + theme(legend.position = "none"),
    cor_plot + theme(axis.text = element_text(size = 7), axis.text.x = element_text(angle = 45, vjust = .6), legend.position = "none"),
    ncol = 1, labels = c("B", "C")
  ),
  nrow = 1, labels = "A"
)


SI_fig
#ggsave(SI_fig, filename = paste0(plot_path, "main_figures/SI_figure_EnsembleMHC_benchmarking.pdf"), height = 6, width = 15)
