data_path <- "~/Covid-19/EnsembleMHC-Covid19/datasets/"
Ensemble_PATH <- "~/Covid-19/EnsembleMHC-Covid19"

library(ggthemes)
library(ggpubr)
library(ggplot2)
library(data.table)
library(stringr)
min_norm <- function(x) {
  (max(x) - x) / (max(x) - min(x))
}
std <- function(x) {
  (x - mean(x)) / sd(x)
}
library(latex2exp)
library(lubridate)
library(tidyr)
library(viridis)
library(ggrepel)
library(reshape2)
library(pwr)
library(wpp2019)
library(dplyr)
library(parallel)
library(wpp2019)
# load population data
data("pop")
# set algorithm names
algos <- c("mhcflurry_affinity_percentile", "mhcflurry_presentation_score", "MixMHCpred", "netMHC_affinity", "netMHCpan_EL_affinity", "netstab_affinity", "pickpocket_affinity")

# Functions

#
get_day_data <- function(day, death_threshold, MIN_countries) {
  # filter at a set number of deaths
  # set days for time series analysis
  # add days since it reached that point
  day_data <- CoV_data_sel_countries %>%
    slice(which(Deaths >= death_threshold)) %>%
    group_by(country) %>%
    mutate(days = as.numeric((today() - day_correction) - mdy(date))) %>%
    mutate(days = max(days) - days)


  # set the threshold of days so you know when to stop iterating once the number of countries decreases below threshold
  thres <- table(day_data$days)[which(table(day_data$days) >= MIN_countries)]
  thres <- as.numeric(names(thres)[length(thres)])

  # mere the day data with the population data
  score_and_death_pop <- day_data %>% merge(select_population)

  # merge EMP scores and calculate deaths per million
  score_and_death_pop <- score_and_death_pop %>%
    group_by(country, days) %>%
    mutate(pop = `2020` / 1000) %>%
    mutate(death_per_pop = Deaths / pop) %>%
    merge(countryEMPscore) %>%
    select(country, death_per_pop, Deaths, days, All.proteins, Structural.proteins) %>%
    slice(which(Deaths >= death_threshold))

  # return day specific data
  score_and_death_pop %>% slice(which(days == day))
}


get_stats <- function(death_threshold, MIN_countries) {
  # filter at a set number of deaths
  # set days for time series analysis
  # add days since it reached that point
  day_data <- CoV_data_sel_countries %>%
    slice(which(Deaths >= death_threshold)) %>%
    group_by(country) %>%
    mutate(days = as.numeric((today() - day_correction) - mdy(date))) %>%
    mutate(days = max(days) - days)


  # set the threshold of days so you know when to stop iterating once the number of countries decreases below threshold
  thres <- table(day_data$days)[which(table(day_data$days) >= MIN_countries)]
  thres <- as.numeric(names(thres)[length(thres)])

  # mere the day data with the population data
  score_and_death_pop <- day_data %>% merge(select_population)

  # merge EMP scores and calculate deaths per million
  score_and_death_pop <- score_and_death_pop %>%
    group_by(country, days) %>%
    mutate(pop = `2020` / 1000) %>%
    mutate(death_per_pop = Deaths / pop) %>%
    merge(countryEMPscore) %>%
    select(country, death_per_pop, Deaths, days, All.proteins, Structural.proteins) %>%
    slice(which(Deaths >= death_threshold))

  cor_data <- lapply(1:thres, function(i) {
    # create subset of data at day i
    print(i)
    tmp <- score_and_death_pop %>% slice(which(days == i))
    list(spearmanConfInt(tmp$Structural.proteins, tmp$death_per_pop),
      spearmanConfInt(tmp$All.proteins, tmp$death_per_pop),
      num_of_country = nrow(tmp)
    )
  })

  data.frame(rbind(
    data.frame(prot_group = "structural", t(sapply(cor_data, function(s) s[[1]])), day = 1:thres),
    data.frame(prot_group = "full", t(sapply(cor_data, function(s) s[[2]])), day = 1:thres)
  ), death_threshold = death_threshold, Countries = sapply(cor_data, function(s) s[[3]]))
}

# EnsembleEMP score function
EMP_score <- function(population) {
  # return dataframe
  data.frame(
    # bind dataframes by row
    do.call(
      rbind,
      # iterate through populations
      lapply(population, function(df) {
        # iterate through select alleles in that population
        # note: because of the 95% allele threshold, there will be 50-52 alleles
        # this is accounted for by dividing by the number of alleles
        apply(t(sapply(as.character(str_remove(df$HLA, "\\*")), function(i) {

          # calculate the EMP score based on all proteins or on structural proteins
          sapply(c("All proteins", "Structural proteins"), function(j) {
            # pf is the peptide fraction
            pf <- HLA_count$pep_frac[which(HLA_count$HLA == i & HLA_count$source == j)]
            # score is the peptide fraction weighted by the normalized expression of that allele
            score <- (df$pop_freq[which(df$HLA == i)] * pf)
            # if the co
            # if(length(score)==0){0}else{score}
          })
        })), 2, sum) / length(df$HLA)
      })
    ),
    stringsAsFactors = F
  )
}


# spearman CI Jackknife  Euclidean likelihood modified from spearmanCI pacakge
#' article{de2012jackknife,
#'   title={Jackknife Euclidean Likelihood-Based Inference for Spearman's Rho},
#'     author={de Carvalho, Miguel and Marques, Filipe J},
#'     journal={North American Actuarial Journal},
#'     volume={16},
#'     number={4},
#'     pages={487--492},
#'     year={2012},
#'     publisher={Taylor \& Francis}
#'   }
spearmanConfInt <- function(x, y, level = 0.95) {
  nx <- length(x)
  ny <- length(y)
  if (is.vector(x) == FALSE | is.vector(x) == FALSE | nx !=
    ny) {
    stop("x and y must be vectors of the same length")
  }
  if (sum(is.na(x)) != 0 | sum(is.na(y)) != 0) {
    stop("missing values are not allowed")
  }
  spear <- cor.test(x, y, method = "spearman")
  U <- spear$estimate
  P <- spear$p.value
  n <- nx
  Z <- as.double()
  for (i in 1:n) Z[i] <- n * U - (n - 1) * cor(x[-i], y[-i], method = "spearman")

  s <- function(theta) 1 / n * sum((Z - theta)^2)
  g <- function(theta) n * (U - theta)^2 / s(theta) - qchisq(level, 1)
  l <- uniroot(g, interval = c(-1, U), tol = 0.1 * 10^{
    -10
  }, extendInt = "yes")$root
  if (l < -1) {
    l <- -1
  }
  # due to the use of the bisection method, sometime teh interval needs to be extended slightly in order to find a root
  # the result of this is lower bound correlation of less than -1
  # stack overflow https://stackoverflow.com/questions/38961221/uniroot-solution-in-r
  u <- uniroot(g, interval = c(U, 1), tol = 0.1 * 10^{
    -10
  })$root
  c(estimate = U, lower = l, upper = u, p_value = P)
}


# this is taken from the JHU github page: https://github.com/CSSEGISandData/COVID-19
# load JHU data
all_data_clean <- read.csv(paste0(data_path, "JHU_data_clean.csv"))


# day correction for when the data ends.
day_correction <- as.numeric(today() - mdy("4/10/2020"))

# read in the matrix of all proteins and structural proteins that were calculated previously
df_all_proteins <- read.csv(paste0(data_path, "identified_peptides_all_proteins_summarise_HLA_protein_counts.csv"))
df_struct <- read.csv(paste0(data_path, "identified_peptides_all_structural_proteins_only_summarise_HLA_protein_counts.csv"))

# creates the count matrix for peptide scoring. This is the total number of peptides predicted for each allele with respect to protein set
# the peptide frraction is also calculated
HLA_count <- rbind(df_all_proteins %>% mutate(source = "All proteins") %>%
  select(HLA, count = total_count_all_proteins, source), df_struct %>% mutate(source = "Structural proteins") %>%
  select(HLA, count = total_count_struct, source)) %>%
  group_by(source) %>%
  mutate(pep_frac = count / sum(count))



all_data_clean <- read.csv(paste0(data_path, "JHU_data_clean.csv"), stringsAsFactors = F)

all_HLA_data <- read.csv(paste0(data_path, "HLA_data_manuscript.csv"), stringsAsFactors = F, header = T)
# countries with at least 1 case reported case by 04/10
confirmed_case_threshold <- 1
test_set <- unique(all_data_clean$Country.Region[which(all_data_clean$Confirmed >= confirmed_case_threshold)])
# Duplicate variable
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
    summarise(allele_count = sum(eff_pop)) %>% # group aggreate by allele
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



# apply the peptide score function
countryEMPscore <- EMP_score(allele_freq_list)
# add a column for the names of the countries
countryEMPscore$country <- row.names(countryEMPscore)

# fix the column name so you can merge later
colnames(CoV_data_sel_countries)[1] <- "country"




# make dataframe of 2020 populations by country
select_population <- pop %>% select(name, "2020")

# fix the offending names so that they match with the population data
colnames(select_population)[1] <- "country"
select_population$country <- as.character(select_population$country)
select_population$country[str_which(select_population$country, "Korea")][2] <- "South Korea"
select_population$country[str_which(select_population$country, "Iran")] <- "Iran"
select_population$country[str_which(select_population$country, "Russia")] <- "Russia"
select_population$country[str_which(select_population$country, "United States of ")] <- "US"
select_population$country[str_which(select_population$country, "United Kingdom")] <- "UK"
select_population$country[str_which(select_population$country, "Taiw")] <- "Taiwan"
select_population$country[str_which(select_population$country, "Hong")] <- "Hong Kong"


all_stats <- mclapply(1:100, function(death_thres) {
  get_stats(death_thres, 8)
}, mc.cores = 10)
all_s <- do.call(rbind, all_stats)

CI_plot <- all_s %>%
  slice(which(prot_group == "structural")) %>%
  ggplot(aes(x = day, y = estimate.rho, color = prot_group)) +
  geom_line(color = "#00ACED") +
  geom_ribbon(aes(ymax = upper, ymin = lower), alpha = .2) +
  theme_classic() +
  facet_wrap(death_threshold ~ .) +
  theme(legend.position = "none") +
  ggtitle("SARS-CoV-2 structural proteins") +
  xlab("days from death threshold")

ggsave(CI_plot, file = paste0(Ensemble_PATH, "/plots/SI_figures/SI_SARS-CoV-2_structural_CI.pdf"), width = 11, height = 11)

CI_plot <- all_s %>%
  slice(which(prot_group == "full")) %>%
  ggplot(aes(x = day, y = estimate.rho, color = prot_group)) +
  geom_line(color = "#388659") +
  geom_ribbon(aes(ymax = upper, ymin = lower), alpha = .2) +
  theme_classic() +
  facet_wrap(death_threshold ~ .) +
  theme(legend.position = "none") +
  ggtitle("Full SARS-CoV-2 proteome") +
  xlab("days from death threshold")

ggsave(CI_plot, file = paste0(Ensemble_PATH, "/plots/SI_figures/SI_Full_SARS-CoV-2_CI.pdf"), width = 11, height = 11)


country_count <- all_s %>% ggplot(aes(x = day, y = Countries)) +
  geom_line() +
  theme_classic() +
  facet_wrap(death_threshold ~ .) +
  xlab("days from death threshold") +
  ylab("number of countries") +
  ggtitle("number of countries in each correlation as a function of time")

ggsave(country_count, file = "SI_country_count_by_day.pdf", width = 11, height = 11)

write.csv(all_s, file = "SI_table_2_cor_data.csv", row.names = F)