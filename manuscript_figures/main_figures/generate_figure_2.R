# read in paths
source("~/Covid-19/EnsembleMHC-Covid19/manuscript_figures/set_paths.R")

library(ggthemes)
library(ggpubr)
library(ggplot2)
library(dplyr)
library(parallel)
library(data.table)
library(stringr)
library(GGally)
library(wesanderson)
min_norm <- function(x) {
  (max(x) - x) / (max(x) - min(x))
}
std <- function(x) {
  (x - mean(x)) / sd(x)
}
library(latex2exp)
library(ggridges)
library(reshape2)


# set algorithm name vector
algos <- c("mhcflurry_affinity_percentile", "mhcflurry_presentation_score", "MixMHCpred", "netMHC_affinity", "netMHCpan_EL_affinity", "netstab_affinity", "pickpocket_affinity")

# load parameterization summary matrix. This is the stored algorithm and allele specific score and FDR
load(paste0(data_path, "P_sum_median_1000_boot.R"))

# load colors for the change plot
c25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)

# read in the coronavirus predictions
corona_data <- fread(pred_peptides, stringsAsFactors = F)

# calculate peptide probabilities
# start by iterating through every unique HLA
pep_probs <- lapply(unique(corona_data$HLA), function(w) {
  # create a tmp variable consisting of the corona virus predictions for one allele
  tmp <- corona_data %>%
    dplyr::slice(which(corona_data$HLA == w)) %>%
    data.frame()
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
  # combine all of the prob_list elements in one dataframe, calculate the products of the FDRs associated with each detecting algorithms
  # merge with information regarding the gene and HLA
  prob_combo <- do.call(rbind, prob_list) %>%
    data.frame() %>%
    group_by(peptide) %>%
    summarise(prob = prod(prob)) %>%
    merge(tmp[, c("peptide", "gene", "HLA")])
})



# bind all elements of  pep_prob list generated in the previous step
all_counts <- do.call(rbind, pep_probs)



# make vector of all unique alleles
sel_alleles <- unique(all_counts$HLA)

# probabiliyty threshold
prob_threshold <- .05

# Summarize by allele and filter for peptides that are less than or equal to 5%
# return a count of number of peptides assigned to each protein with respect to each allele 
df_all_proteins <- all_counts %>%
  slice(which(prob <= prob_threshold)) %>%
  group_by(HLA, gene) %>%
  summarise(count = length(gene))


# convert gene name to character so orf1ab name can be corrected
df_all_proteins$gene <- as.character(df_all_proteins$gene)
# rename orf1ab  gene to make it match other gene name formats
df_all_proteins$gene[which(df_all_proteins$gene == "orf1ab")] <- "ORF1ab"
# convert back to a factor
df_all_proteins$gene <- factor(df_all_proteins$gene)
# add name protein group name to the matrix
df_all_proteins$type <- "All proteins"
# make bar graph based on the all proteins. This is figure 2A
p1 <- df_all_proteins %>% ggplot(aes(x = factor(HLA, levels = rev(unique(HLA))), y = count, fill = gene)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  ggtitle("full SARS-CoV-2 protome") +
  scale_fill_manual(values = c25[1:10]) +
  theme_linedraw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("number of peptides") +
  xlab("HLA")

# create another bar plot of only structural proteins
df_struct <- all_counts %>%
  slice(which(prob <= prob_threshold)) %>%
  slice(which(gene %in% c("E", "N", "M", "S"))) %>%
  group_by(HLA, gene) %>%
  summarise(count = length(gene))
# set name of protein group
df_struct$type <- "Structural proteins"
# for plotting, HLAs with no predicted structural peptides are set to 0. This is to avoid unessessary dropping of HLAs
# This is done by looking for missing alleles, assigning a value of zero to an arbitrary structural protein, in this case the E protein
# these additional pseudo counts are then combined with the df_struct matrix
df_struct <- rbind(df_struct %>% data.frame(), data.frame(HLA = unique(df_all_proteins$HLA)[-which(unique(df_all_proteins$HLA) %in% unique(df_struct$HLA))], gene = "E", count = 0, type = "Structural protein"))
# factor the HLA names
df_struct$HLA <- factor(df_struct$HLA, levels = rev(unique(df_struct$HLA)[order(unique(df_struct$HLA))]))
# plot the bar graph for peptide-allele distrubtion for only structural proteins. Figure 2B
p2 <- df_struct %>% ggplot(aes(x = HLA, y = count, fill = gene)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  ggtitle("SARS-CoV-2 structural proteins") +
  scale_fill_manual(values = c25[c(1:3, 10)]) +
  theme_linedraw() +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("number of peptides") +
  scale_y_continuous(breaks = seq(0, 12, 2))



# calculate the standev of each allele with respect to all proteins
df_all_w_std <- df_all_proteins %>%
  group_by(HLA) %>%
  summarise(total_count_all_proteins = sum(count)) %>%
  mutate(std_all = std(total_count_all_proteins)) %>%
  mutate(pep_frac_all = total_count_all_proteins / sum(total_count_all_proteins))

#print some stats
median(df_all_w_std$total_count_all_proteins)
max(df_all_w_std$total_count_all_proteins)
min(df_all_w_std$total_count_all_proteins)
IQR(df_all_w_std$total_count_all_proteins)


# combine all structural protein predictions and calculate standard deviation
df_struct_w_std <- df_struct %>%
  group_by(HLA) %>%
  summarise(total_count_struct = sum(count)) %>%
  mutate(std_struct = std(total_count_struct)) %>%
  mutate(pep_frac_struct = total_count_struct / sum(total_count_struct))

#print some stats
summary(df_struct_w_std$total_count_struct)
IQR(df_struct_w_std$total_count_struct)


# merge both sets
all_combo <- df_all_w_std %>% merge(df_struct_w_std)

# The absolute difference in peptide fraction by allele
all_combo$delta <- all_combo$pep_frac_all - all_combo$pep_frac_struct

# median peptide fraction
med_pep_frac <- median(c(all_combo$pep_frac_all, all_combo$pep_frac_struct))
# covert all alleles that change by less than the median peptide fraction to NC. NC = no significant change
all_combo$HLA[which(abs(all_combo$delta) < med_pep_frac)] <- "NC"
# set alpha levels by allele based on change. if they change by  > median peptide fraction they have alpha of one, else alpha of .2
all_combo$alpha <- 1
all_combo$alpha[which(all_combo$HLA == "NC")] <- .2
all_combo$HLA <- factor(all_combo$HLA)
# rename change columns
colnames(all_combo)[c(4, 7)] <- c("full SARS-CoV-2 proteome", "SARS-CoV-2 structural proteins")

# make parallel plot
p3 <- ggparcoord(all_combo, columns = c(4, 7), groupColumn = "HLA", scale = "globalminmax", alphaLines = "alpha") + scale_color_manual(values = c(c25[c(1:5, 7)], "black")) +
  scale_alpha(guide = "none") + ylab("peptide fraction") +
  xlab("") + theme_linedraw() + theme(legend.key.size = unit(2, "mm")) + ggtitle("relative change in peptide fraction")


garg_1 <- ggarrange(p1, p2, common.legend = T, legend = "bottom")

ggarrange(garg_1, p3)

# ggsave(ggarrange(garg_1, p3), filename = paste0(Ensemble_PATH, "/plots/main_figures/Figure_2_pep_frac.pdf"), height = 12, width = 15)

# the following files will be useful down the line
# this will write the file to the dataset directory. It will be used in subsequnt scripts

# all peptides before FDR filter
write.csv(all_counts, file = paste0(Ensemble_PATH, "/datasets/all_peptides_prefilter.csv"), row.names = F)

# all peptides after FDR filter
write.csv(df_all_proteins, file = paste0(Ensemble_PATH, "/datasets/all_peptides_passing_score_filter_summarized_by_FDR.csv"), row.names = F)
# the dataframe of total number of peptide per allele for all proteins
write.csv(df_all_w_std, paste0(data_path, "/identified_peptides_all_proteins_summarise_HLA_protein_counts.csv"), row.names = F)
# teh dataframe of total number of peptide per allele for structural proteins
write.csv(df_struct_w_std, paste0(data_path, "/identified_peptides_all_structural_proteins_only_summarise_HLA_protein_counts.csv"), row.names = F)

# these are the same as above but retain the protein level information
write.csv(df_all_proteins, paste0(data_path, "/identified_peptides_all_proteins_summarise_w_gene_counts.csv"), row.names = F)
write.csv(df_struct, paste0(data_path, "/identified_peptides_all_structural_proteins_only_summarise_HLA_protein_w_gene_counts.csv"), row.names = F)
