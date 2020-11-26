source("~/Covid-19/EnsembleMHC-Covid19/manuscript_figures/set_paths.R")
library(ggpubr)
library(ggplot2)
library(dplyr)
library(data.table)
library(stringr)
library(patchwork)
algos <- c("mhcflurry_affinity_percentile", "mhcflurry_presentation_score", "MixMHCpred", "netMHC_affinity", "netMHCpan_EL_affinity", "netstab_affinity", "pickpocket_affinity")


all_proteins <- fread(paste0(data_path, "all_peptides_prefilter.csv"))
p1 <- all_proteins %>% ggplot(aes(x = prob)) +
  geom_density() +
  geom_vline(xintercept = .05, color = "red") +
  theme_pubclean() +
  xlab("peptide FDR")

all_proteins <- all_proteins %>% slice(which(prob <= .05))

struct_tab <- table(nchar(all_proteins$peptide[which(all_proteins$gene %in% c("E", "N", "M", "S"))]))
all_tab <- table(nchar(all_proteins$peptide))

p2 <- data.frame(len = as.numeric(names(struct_tab)), count = as.numeric(struct_tab) / sum(struct_tab)) %>%
  ggplot(aes(len, count)) +
  geom_bar(stat = "identity", fill = "#00ACED") +
  xlab("peptide length") +
  ylab("% of peptides") +
  ggtitle("SARS-CoV-2 structural proteins") +
  theme_pubclean() +
  scale_x_continuous(breaks = 8:14)

p3 <- data.frame(len = as.numeric(names(all_tab)), count = as.numeric(all_tab) / sum(all_tab)) %>%
  ggplot(aes(len, count)) +
  geom_bar(stat = "identity", fill = "#ED0036") +
  xlab("peptide length") +
  ylab("% of peptides") +
  ggtitle("full SARS-CoV-2 proteome") +
  theme_pubclean() +
  scale_x_continuous(breaks = 8:14)


garg <- p1 + {
  p2 / p3
}
garg
#ggsave(garg, filename = paste0(Ensemble_PATH, "plots/SI_figures/SI_peplen_prod_dist.pdf"), width = 12, height = 8)
