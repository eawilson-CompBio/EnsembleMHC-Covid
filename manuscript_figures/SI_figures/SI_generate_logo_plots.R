Ensemble_PATH <- "~/Covid-19/EnsembleMHC-Covid19"
data_path <- "~/Covid-19/EnsembleMHC-Covid19/datasets/"
library(Biostrings)
library(ggplot2)
library(dplyr)
library(data.table)
library(stringr)
library(ggseqlogo)


all_counts <- fread(paste0(data_path, "all_peptides_prefilter.csv"))

peps <- all_counts %>%
  slice(which(prob <= .05)) %>%
  group_split(HLA)
# grab peptide, if peptide is logner than 9 remove AA 9:length-1 if peptide is 8 add random AA at 8 move 8 to 9
peps <- lapply(peps, function(P) {
  sapply(P$peptide, function(H) {
    if (nchar(H) > 9) {
      paste(str_split(H, "", simplify = T)[c(1:8, nchar(H))], collapse = "")
    } else if (nchar(H) < 9) {
      tmp <- str_split(H, "", simplify = T)
      paste(paste(tmp[1:7], collapse = ""), AA_STANDARD[sample(1:20, 1)], tmp[8], sep = "")
    } else {
      H
    }
  })
})
names(peps) <- unique(all_counts$HLA)
# remove bottom 95% (>4 peptides)
thres <- quantile(as.numeric(table(all_counts %>% slice(which(prob <= .05)) %>% pull(HLA))), prob = seq(0, 1, .05))[2]
peps <- peps[which(sapply(peps, length) > thres)]

gg <- ggseqlogo(peps, facet = "wrap")


ggsave(gg, filename = paste0(Ensemble_PATH, "/plots/SI_figures/SI_logoplot.pdf"), width = 15, height = 15)
