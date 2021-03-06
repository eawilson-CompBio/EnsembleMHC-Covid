# read in paths
source("~/Covid-19/EnsembleMHC-Covid19/manuscript_figures/set_paths.R")
library(ggplot2)
library(dplyr)
library(parallel)
library(data.table)
library(patchwork)

#-----------------------------------------------------------------------------------------------------
# note: Because the distribution is sampled, it will not directly resemble the plot in the manuscript.
# However, it will still demonstrate the main point of these plots
#-----------------------------------------------------------------------------------------------------
df_all_w_std <- read.csv(paste0(data_path, "/identified_peptides_all_proteins_summarise_HLA_protein_counts.csv"))
df_struct_w_std <- read.csv(paste0(data_path, "/identified_peptides_all_structural_proteins_only_summarise_HLA_protein_counts.csv"))
all_sample <- round(rnorm(n = 52, mean = 16, sd = 1))
ks.test(all_sample, df_all_w_std$total_count_all_proteins)


p1 <- data.frame(HLA = factor(df_all_w_std$HLA, levels = rev(df_all_w_std$HLA)), count = all_sample) %>% ggplot(aes(x = HLA, y = count)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  ylab("number of peptides") +
  ggtitle("full SARS-CoV-2 proteome: uniform distrubtion") +
  theme_classic()

p2 <- data.frame("pep" = all_sample) %>% ggplot(aes(pep)) +
  geom_density() +
  xlab("number of peptides") +
  theme_classic() +
  ggtitle("Full SARS-CoV-2 proteome")

struct_sample <- round(rnorm(n = 52, mean = 2, sd = .5))
p3 <- data.frame(HLA = factor(df_all_w_std$HLA, levels = rev(df_all_w_std$HLA)), count = struct_sample) %>% ggplot(aes(x = HLA, y = count)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  ylab("number of peptides") +
  ggtitle("SARS-CoV-2 structural proteins: uniform distrubtion") +
  theme_classic()


p4 <- data.frame("pep" = struct_sample) %>% ggplot(aes(pep)) +
  geom_density() +
  xlab("number of peptides") +
  scale_x_continuous(breaks = seq(0, max(struct_sample), 1)) +
  theme_classic() +
  ggtitle("SARS-CoV-2 structural proteins")


#SI_even_dist <- p1 + p2 + p3 + p4 + plot_layout(ncol = , byrow = F, heights = c(.7, .25))

SI_even_dist <- p2 / p4

#ggsave(filename = paste0(Ensemble_PATH, "/plots/SI_figures/SI_even_distrubtion.pdf"), SI_even_dist, width = 14, height = 10.667)
ggsave(filename = paste0(Ensemble_PATH, "/plots/SI_figures/SI_even_distrubtion.pdf"), SI_even_dist, width = 2.8, height = 2.8)
