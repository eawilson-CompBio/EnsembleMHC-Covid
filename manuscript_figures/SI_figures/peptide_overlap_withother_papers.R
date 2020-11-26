# read in paths
source("~/Covid-19/EnsembleMHC-Covid19/manuscript_figures/set_paths.R")

library(patchwork)
library(data.table)
library(pracma)
library(parallel)
library(dplyr)
library(stringr)
library(ggplot2)

all_corona_peps <- fread(paste0(data_path,"/predicted_corona_peptides.csv")) %>%
  data.frame() %>%
  mutate(pickpocket_affinity = 50000^(1 - pickpocket_affinity)) %>%
  slice(which(gene%in%c("S","N","M","E")))

all_corona_peptides <- unique(all_corona_peps$peptide)

Quadeer <- read.csv(paste0(validation_peps,"/summarized_immunogeic_peptides_Quadeer_et_al.csv")) 
Quadeer_peptides <- Quadeer$Epitope
Quadeer_peptides <- Quadeer_peptides[which(Quadeer_peptides%in%all_corona_peptides)]

snyder <- read.csv(paste0(validation_peps,"/immunocode_peptides.csv")) 
snyder_peptides <- unlist(lapply(unique(snyder$Amino.Acids),function(s)str_split(s,pattern = ",",simplify = T)))
snyder_peptides <- snyder_peptides[which(snyder_peptides%in%all_corona_peptides)]


ferretti <- read.csv(paste0(validation_peps,"/Ferretti_et_al.csv"))
ferretti_peptides <- ferretti$Full_Peptide
ferretti_peptides <- ferretti_peptides[which(ferretti_peptides%in%all_corona_peptides)]

nelde <- read.csv(paste0(validation_peps,"/nelde_et_al_peptides.csv"))
nelde <-nelde[-intersect(str_which(nelde$SARS_exp,"^0/"),str_which(nelde$healthy_never_exp,"^0/")),]
nelde_peptides <- str_remove(nelde$peptide,"\xca")
nelde_peptides <- nelde_peptides[which(nelde_peptides%in%all_corona_peptides)]



all_peptides <- list(Snyder_et_al=snyder_peptides,Quadeer_et_al=Quadeer_peptides,Ferretti_et_al=ferretti_peptides,Nelde_et_al=nelde_peptides)


Wilson_peptides <- read.csv(paste0(data_path,"/all_peptides_prefilter.csv")) %>%
  slice(which(prob <= .05)) %>%
  slice(which(gene %in% c("S", "N", "M", "E"))) %>%
  pull(peptide) %>%
  unique()


peptide_sum <- do.call(rbind, lapply(names(all_peptides), function(x) {
  data.frame(table(Wilson_peptides %in% all_peptides[[x]]), pep_set = x)
})) %>% slice(which(Var1==T))

all_valid<-data.frame(table(Wilson_peptides %in% unique(do.call("c", all_peptides))), pep_set = "total validated")


snyder <- read.csv(paste0(validation_peps,"/immunocode_peptides.csv")) 
snyder_peptides <- unique(snyder$Amino.Acids)[which(sapply(unique(snyder$Amino.Acids), function(s) length(str_split(s, pattern = ",", simplify = T))<2))]
snyder_peptides <- snyder_peptides[which(snyder_peptides%in%all_corona_peptides)]
all_peptides <- list(Snyder_et_al=snyder_peptides,Ferretti_et_al=ferretti_peptides,Nelde_et_al=nelde_peptides,Quadeer_peptides)

summary_dat <- rbind(
  all_valid,
  data.frame(table(Wilson_peptides %in% unique(do.call("c", all_peptides))), pep_set = "total validated (no pools)")
) %>% slice(which(Var1==T))

ggthemr::ggthemr("fresh")
p1 <- peptide_sum %>% 
  ggplot(aes(x = pep_set, y = Freq)) +
  geom_bar(stat = "identity", fill = "#00ACED") +
  ggtitle("validated peptides by study") +
  theme(legend.position = "top") +
  labs(fill = "observed in 108 structural protein peptide set data set") +
  geom_text(data = peptide_sum, mapping = aes(label = Freq, y = peptide_sum$Freq), color = "black",vjust= -.1)

p2 <- summary_dat %>% ggplot(aes(x = pep_set, y = Freq)) +
  geom_bar(stat = "identity", fill = "#ED4100") +
  ggtitle("peptide summary") +
  theme(legend.position = "top") +
  labs(fill = "observed in 108 structural protein peptide set data set") +
  geom_text(data = summary_dat , aes(label = Freq, y = summary_dat$Freq), color = "black",vjust= -.1)


p1 + p2 + plot_layout(guides = "collect") & theme(legend.position = "top") & ylab("number of peptides") & xlab("")

