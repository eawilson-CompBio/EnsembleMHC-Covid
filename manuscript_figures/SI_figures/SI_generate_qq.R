Ensemble_PATH <- "~/Covid-19/EnsembleMHC-Covid19"
dataset_path <- paste0(Ensemble_PATH, "/datasets/")
library(ggpubr)
library(dplyr)
library(data.table)
library(stringr)
library(lubridate)
library(patchwork)
# this is taken from the JHU github page: https://github.com/CSSEGISandData/COVID-19
# load JHU data 1/22 - 4/10
all_data_clean <- fread(paste0(dataset_path, "JHU_data_clean.csv"))


# day correction for when the data ends.
day_correction <- as.numeric(today() - mdy("4/10/2020"))



algos <- c("mhcflurry_affinity_percentile", "mhcflurry_presentation_score", "MixMHCpred", "netMHC_affinity", "netMHCpan_EL_affinity", "netstab_affinity", "pickpocket_affinity")


all_data_clean <- fread(paste0(dataset_path, "JHU_data_clean.csv"))
# day correction
day_correction <- as.numeric(today() - mdy("4/10/2020"))

# read in the matrix of all proteins and structural proteins that were calculated previously
df_all_proteins <- read.csv(paste0(dataset_path, "identified_peptides_all_proteins_summarise_HLA_protein_counts.csv"))
df_struct <- read.csv(paste0(dataset_path, "identified_peptides_all_structural_proteins_only_summarise_HLA_protein_counts.csv"))

# creates the count matrix for peptide scoring. This is the total number of peptides predicted for each allele with respect to protein set
# the peptide frraction is also calculated
HLA_count <- rbind(df_all_proteins %>% mutate(source = "All proteins") %>%
  select(HLA, count = total_count_all_proteins, source), df_struct %>% mutate(source = "Structural proteins") %>%
  select(HLA, count = total_count_struct, source)) %>%
  group_by(source) %>%
  mutate(pep_frac = count / sum(count))


# EnsembleEMP score function
EMP_score <- function(population) {
  # return dataframe
  data.frame(
    # bind dataframes by row
    do.call(
      rbind,
      # iterate through population list
      lapply(population, function(df) {
        # iterate through select alleles in that population
        # note: because of the 95% allele threshold, there will be 50-52 alleles fpr each country
        # this is accounted for by dividing by the number of alleles
        apply(t(sapply(as.character(str_remove(df$HLA, "\\*")), function(i) {

          # calculate the EMP score based on all proteins or on structural proteins
          sapply(c("All proteins", "Structural proteins"), function(j) {
            # pf is the peptide fraction
            pf <- HLA_count$pep_frac[which(HLA_count$HLA == i & HLA_count$source == j)]
            # score is the peptide fraction weighted by the normalized expression of that allele
            score <- (df$pop_freq[which(df$HLA == i)] * pf)
          })
        })), 2, sum) / length(df$HLA)
      })
    ),
    stringsAsFactors = F
  )
}


# read HLA data
all_HLA_data <- fread(paste0(dataset_path, "HLA_data_manuscript.csv"), header = T)
# countries with at least 1 case reported case by 04/10
test_set <- unique(all_data_clean$Country.Region[which(all_data_clean$Deaths >= 1)])
# save country names
country_list <- test_set

# convert to list and fix country names or add alternative names
test_set <- lapply(test_set, function(w) {
  w
})
test_set[[which(country_list == "Czechia")]] <- "Czech"
test_set[[which(country_list == "UK")]] <- c("United Kingdom|England|Wales$")
test_set[[which(country_list == "US")]] <- "USA"


# set allele vector to the 52 selected alleles
sel_HLA <- unique(HLA_count$HLA)

# normalize the allele frequencies with respect to selected alleles
allele_freq_list <- lapply(test_set, function(w) {
  name <- w
  print(name)
  w <- all_HLA_data[str_which(all_HLA_data$population, w), ]
  colnames(w) <- c("HLA", "pop", "allele_freq", "n")
  w$HLA <- str_remove(w$HLA, "\\*")
  total_n <- w %>%
    slice(which(HLA %in% sel_HLA)) %>%
    mutate(country = name) %>%
    select(pop, n) %>%
    unique() %>%
    mutate(n_total = sum(n)) %>%
    pull(n_total) %>%
    unique()


  tmp <- w %>%
    mutate(eff_pop = allele_freq * n * 2) %>% # convert to allele count
    slice(which(HLA %in% sel_HLA)) %>% # isolate the 52 alleles
    group_by(HLA) %>%
    summarise(allele_count = sum(eff_pop)) %>% # group and aggreate allele count by 52 alleles
    mutate(pop_freq = allele_count / sum(allele_count)) %>% # normalize allele count
    mutate(country = name) %>%
    mutate(total_n = total_n) # add country name and total number of HLA typed individual
})

# name list
names(allele_freq_list) <- country_list

# get rid of empty lists
allele_freq_list <- allele_freq_list[which(sapply(allele_freq_list, function(w) {
  nrow(w) > 0
}))]

# filter for populations with at least 1000 HLA typed individuals at 4 digit resolution
allele_freq_list <- allele_freq_list[which(sapply(allele_freq_list, function(w) {
  unique(w$total_n) > 1000
}))]

# filter for number of select alleles
allele_freq_list <- allele_freq_list[which(sapply(allele_freq_list, nrow) / 52 > .95)]
# get coroan virus data that correspondes to the selected alleles
CoV_data_sel_countries <- all_data_clean[which(all_data_clean$Country.Region %in% names(allele_freq_list)), ]

# fix the column name so you can merge later
colnames(CoV_data_sel_countries)[1] <- "country"




all_big <- do.call(rbind, allele_freq_list)
sample <- allele_freq_list[1]
sim_pops <- lapply(1:10000, function(rep) {
  tmp <- sample
  tmp[[1]]$pop_freq <- all_big$pop_freq[sample(1:nrow(all_big), 52)]
  tmp[[1]] <- tmp[[1]] %>% mutate(pop_freq = pop_freq / sum(pop_freq))
})

out <- EMP_score(sim_pops)

p1 <- ggqqplot(out$Structural.proteins) + ylab("sampled EnsembleMHC population score") + ggtitle("sampled SARS-CoV-2 Structural proteins Q-Q plot")
p2 <- ggqqplot(out$All.proteins) + ylab("sampled EnsembleMHC population score") + ggtitle("sampled full SARS-CoV-2 proteome Q-Q plot")

pop_big <- EMP_score(allele_freq_list)
# add a column for the names of the countries
pop_big$country <- row.names(pop_big)


library(wpp2019)

data("pop")


select_population <- pop %>% select(name, "2020")

colnames(select_population)[1] <- "country"
select_population$country <- as.character(select_population$country)
select_population$country[str_which(select_population$country, "Korea")][2] <- "South Korea"
select_population$country[str_which(select_population$country, "Iran")] <- "Iran"
select_population$country[str_which(select_population$country, "Russia")] <- "Russia"
select_population$country[str_which(select_population$country, "United States of ")] <- "US"

day_data <- CoV_data_sel_countries %>%
  slice(which(Deaths >= 1)) %>%
  group_by(country) %>%
  mutate(days = as.numeric((today() - day_correction) - mdy(date))) %>%
  mutate(days = max(days) - days)
# a small number here is just a waste of time at some point


score_and_death_pop <- day_data %>% merge(select_population)

score_and_death_pop <- score_and_death_pop %>%
  group_by(country, days) %>%
  mutate(popM = `2020` / 1000) %>%
  mutate(death_per_pop = Deaths / popM) %>%
  merge(pop_big) %>%
  select(country, death_per_pop, Deaths, days, All.proteins, Structural.proteins)


p3 <- ggqqplot(score_and_death_pop$death_per_pop) + ylab("Death rate per million") + ggtitle("Deaths per million Q-Q plot")

SI_AB <- p1 + p2
SI_C <- p3

ggsave(SI_AB, paste0(Ensemble_PATH, "/plots/SI_figures/SI_QQ_AB.pdf"))
ggsave(SI_C, paste0(Ensemble_PATH, "/plots/SI_figures/SI_QQ_C.pdf"))