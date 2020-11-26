# read in paths
source("~/Covid-19/EnsembleMHC-Covid19/manuscript_figures/set_paths.R")

library(ggplot2)
library(data.table)
library(stringr)
library(lubridate)
library(tidyr)
library(viridis)
library(reshape2)
library(wpp2019)
library(dplyr)
library(parallel)
library(ggthemr)
library(patchwork)
library(LaplacesDemon)
ggthemr("fresh")

##functions


#EnsembleEMP score function
EMP_score <- function(population,HLA_count) {
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
          })
        })), 2, sum) / length(df$HLA)
      })
    ),
    stringsAsFactors = F
  )
}


#calculate correlation as a function of time
death_threshold_specific_corr <- function(death_threshold, MIN_countries,countryEMPscore) {
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
    tmp <- score_and_death_pop %>% slice(which(days == i))
    list(
      cor_struct = cor(tmp$Structural.proteins, tmp$death_per_pop, method = "spearman"),
      cor_all_prot = cor(tmp$All.proteins, tmp$death_per_pop, method = "spearman"),
      tab = tmp,
      day = i,
      p_struct = cor.test(tmp$Structural.proteins, tmp$death_per_pop, method = "spearman")[[3]],
      p_all_prot = cor.test(tmp$All.proteins, tmp$death_per_pop, method = "spearman")[[3]]
    )
  })

  p_val_struct <- sapply(1:thres, function(w) {
    cor_data[[w]]$p_struct
  })
  p_val_all_prot <- sapply(1:thres, function(w) {
    cor_data[[w]]$p_all_prot
  })
  # make a matrix of protein group specific correlations
  df_corr <- data.frame(days = 1:thres, correlation_all_prot = sapply(cor_data, function(w) {
    w$cor_all_prot
  }), correlation_struct = sapply(cor_data, function(w) {
    w$cor_struct
  }))

  # create the p value matrix
  pvl_mat <- rbind(
    data.frame(days = df_corr$days, source = "Structural_proteins", p_value = p_val_struct),
    data.frame(days = df_corr$days, source = "All_proteins", p_value = p_val_all_prot)
  )

  # renames columns and melt matrix
  colnames(df_corr)[2:3] <- c("All_proteins", "Structural_proteins")
  melted_df <- melt(df_corr, id.vars = c("days"))
  colnames(melted_df)[2:3] <- c("source", "correlation")

  # merge with p value matrix and identify significant correlations
  melted_df <- melted_df %>% merge(pvl_mat)
  melted_df$sig <- 0
  melted_df$sig[which(melted_df$p_value <= .05)] <- 1

  # find the number of countries at each time point regradless of the minimum country threshold
  tab_melt <- table(day_data$days)
  # merge so only considered days meeting the min country remain
  melted_df <- melted_df %>% merge(data.frame(days = as.numeric(names(tab_melt)), num_country = as.numeric(tab_melt)))


  melted_df$confirmed <- death_threshold
  # return matrix
  melted_df
}


# load population data
data("pop")
# set algorithm names
algos <- c("mhcflurry_affinity_percentile", "mhcflurry_presentation_score", "MixMHCpred", "netMHC_affinity", "netMHCpan_EL_affinity", "netstab_affinity", "pickpocket_affinity")

# this is taken from the JHU github page: https://github.com/CSSEGISandData/COVID-19
# load JHU data 1/22 - 4/10
all_data_clean <- fread(paste0(data_path, "JHU_data_clean.csv"))


# day correction for when the data ends.
day_correction <- as.numeric(today() - mdy("4/10/2020"))


# get all peptides predicted by ensembleMHC
df_all_peptides <- read.csv(paste0(data_path, "all_peptides_prefilter.csv"))

#load the HLA list
load(paste0(data_path,"selected_alleles.R"))


# generate the correlations reported in the paper -----------------------

# read in the matrix of all proteins and structural proteins that were calculated previously
df_all_proteins <- read.csv(paste0(data_path, "identified_peptides_all_proteins_summarise_HLA_protein_counts.csv"))
df_struct <- read.csv(paste0(data_path, "identified_peptides_all_structural_proteins_only_summarise_HLA_protein_counts.csv"))

# # creates the count matrix for peptide scoring. This is the total number of peptides predicted for each allele with respect to protein set
# # the peptide frraction is also calculated
 HLA_count <- rbind(df_all_proteins %>% mutate(source = "All proteins") %>%
   select(HLA, count = total_count_all_proteins, source), df_struct %>% mutate(source = "Structural proteins") %>%
   select(HLA, count = total_count_struct, source)) %>%
   group_by(source) %>%
   mutate(pep_frac = count / sum(count))


 real_HLA_count <- HLA_count
 
# read HLA data
all_HLA_data <- fread(paste0(data_path, "HLA_data_manuscript.csv"), header = T)
# countries with at least 1 case reported case by 04/10
test_set <- unique(all_data_clean$Country.Region[which(all_data_clean$Deaths >= 1)])
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
load(paste0(data_path,"selected_alleles.R"))

# normalize the allele frequencies with respect to selected alleles
allele_freq_list <- lapply(test_set, function(w) {
  name <- w
  print(name)
  w <- all_HLA_data[str_which(all_HLA_data$population, w), ]
  colnames(w) <- c("HLA", "pop", "allele_freq", "n")
  w$HLA <- str_remove(w$HLA, "\\*")
  total_n <- w %>%
    slice(which(HLA %in% sel_alleles)) %>%
    mutate(country = name) %>%
    select(pop, n) %>%
    unique() %>%
    mutate(n_total = sum(n)) %>%
    pull(n_total) %>%
    unique()
  
  
  tmp <- w %>%
    mutate(eff_pop = allele_freq * n * 2) %>% # convert to allele count
    slice(which(HLA %in% sel_alleles)) %>% # isolate the 52 alleles
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
# when randomly sampling there is not 52 alleles so This will need to be adjuested accordingly
allele_freq_list <- allele_freq_list[which(sapply(allele_freq_list, nrow) / length(sel_alleles) > .95)]
# allele_freq_list <- allele_freq_list[which(sapply(allele_freq_list, nrow) / 52 > .95)]
# get coroan virus data that correspondes to the selected alleles
CoV_data_sel_countries <- all_data_clean[which(all_data_clean$Country.Region %in% names(allele_freq_list)), ]



# apply the peptide score function
countryEMPscore <- EMP_score(allele_freq_list,HLA_count = HLA_count)
# add a column for the names of the countries
countryEMPscore$country <- row.names(countryEMPscore)

# fix the column name so you can merge later
colnames(CoV_data_sel_countries)[1] <- "country"




# make dataframe of 2020 populations by country
select_population <- pop %>% select(name, "2020")

# some countries names have names that differ than the ones used in the JHU data
# This step changes names that differ to the names used in the JHU data
colnames(select_population)[1] <- "country"
select_population$country <- as.character(select_population$country)
select_population$country[str_which(select_population$country, "Korea")][2] <- "South Korea"
select_population$country[str_which(select_population$country, "Iran")] <- "Iran"
select_population$country[str_which(select_population$country, "Russia")] <- "Russia"
select_population$country[str_which(select_population$country, "United States of ")] <- "US"
select_population$country[str_which(select_population$country, "United Kingdom")] <- "UK"
select_population$country[str_which(select_population$country, "Taiw")] <- "Taiwan"
select_population$country[str_which(select_population$country, "Hong")] <- "Hong Kong"



death_threshold_iter <- mclapply(1:100, function(death_thres) {
  death_threshold_specific_corr(death_thres, 8,countryEMPscore = countryEMPscore)
}, mc.cores = 10)


# create data frame from death threshold correlations
# combine results into a matrix
df_Death <- do.call(rbind, death_threshold_iter)
# normalize days for visualization purposes
norm_df_Death <- do.call(rbind, lapply(death_threshold_iter, function(d) {
  d$days <- ((d$days) - min(d$days)) / ((max(d$days) - min(d$days)))
  d
}))

# create matrix of only significant correlations
# rename corr types
norm_df_Death$source <- as.character(norm_df_Death$source)
norm_df_Death <- norm_df_Death[which(norm_df_Death$source == "Structural_proteins"), ]
norm_df_Death$source[which(norm_df_Death$source == "Structural_proteins")] <- "SARS-CoV-2 structural proteins"
# norm_df_Death$source[which(norm_df_Death$source=="All_proteins")]<-"entire SARS-CoV-2 proteome"
norm_df_Death$source <- as.factor(norm_df_Death$source)
sig <- subset(norm_df_Death, sig == 1)

table(norm_df_Death$sig) / sum(table(norm_df_Death$sig))

values_from_paper <- df_Death %>% select(source,p_value,correlation)
df_death_comp <- df_Death




#subsample data--------------------------------------

# read HLA data
all_HLA_data <- fread(paste0(data_path, "HLA_data_manuscript.csv"), header = T)
# countries with at least 1 case reported case by 04/10
test_set <- unique(all_data_clean$Country.Region[which(all_data_clean$Deaths >= 1)])
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
load(paste0(data_path, "selected_alleles.R"))

# normalize the allele frequencies with respect to selected alleles
allele_freq_list <- lapply(test_set, function(w) {
  name <- w
  print(name)
  w <- all_HLA_data[str_which(all_HLA_data$population, w), ]
  colnames(w) <- c("HLA", "pop", "allele_freq", "n")
  w$HLA <- str_remove(w$HLA, "\\*")
  total_n <- w %>%
    slice(which(HLA %in% sel_alleles)) %>%
    mutate(country = name) %>%
    select(pop, n) %>%
    unique() %>%
    mutate(n_total = sum(n)) %>%
    pull(n_total) %>%
    unique()


  tmp <- w %>%
    mutate(eff_pop = allele_freq * n * 2) %>% # convert to allele count
    slice(which(HLA %in% sel_alleles)) %>% # isolate the 52 alleles
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
# when randomly sampling there is not 52 alleles so This will need to be adjuested accordingly
allele_freq_list <- allele_freq_list[which(sapply(allele_freq_list, nrow) / length(sel_alleles) > .95)]

# get coroan virus data that correspondes to the selected alleles
CoV_data_sel_countries <- all_data_clean[which(all_data_clean$Country.Region %in% names(allele_freq_list)), ]


resample_dist<-list()

for (rev_2 in 1:1000) {

  # Here is subsample the number of peptides from all idenitified peptides to match the number of structural peptides
  # This is to demonstrate that the results are different between the two groups reported in the paper are robust and not due to down sampling
  # one note is that I am duplicating the same sample for all proteins and strucutral proteins groups
  # this is due to a lot of the legacy scripts requiring both of these inputs
  # I want to leave as much intact as possible to avoid any potential changes in analysis that might effect the analysis
  # you can run diff on this script and generate_figure_3.R to determine the differences
  # furthermore, you can check that both all proteins and strucutral proteins proudce the same results

  # get 108 unique randmoly sampled peptides
  peptide_sample <- sample(df_all_peptides %>% slice(which(prob <= .05))  %>% pull(peptide) %>% unique(), 108, replace = F)

  tmp_all_prots <- df_all_peptides %>%
    slice(which(prob <= .05)) %>%
    slice(which(peptide %in% peptide_sample)) %>%
    group_by(HLA) %>%
    summarise(count = length(HLA)) %>%
    mutate(source = "All proteins") %>%
    mutate(pep_frac = count / sum(count))

  if (nrow(tmp_all_prots) < 52) {
    tmp_all_prots <- rbind(tmp_all_prots, data.frame(HLA = sel_alleles[-which(sel_alleles %in% tmp_all_prots$HLA)], count = 0, source = "All proteins", pep_frac = 0))
  }

  tmp_struct_prots <- tmp_all_prots
  tmp_struct_prots$source <- "Structural proteins"

  HLA_count_sample <- rbind(tmp_all_prots, tmp_struct_prots)

  # apply the peptide score function
  countryEMPscore_sample <- EMP_score(allele_freq_list, HLA_count = HLA_count_sample) %>% data.frame()
  # add a column for the names of the countries
  countryEMPscore_sample$country <- row.names(countryEMPscore_sample)

  # fix the column name so you can merge later
  colnames(CoV_data_sel_countries)[1] <- "country"




  # make dataframe of 2020 populations by country
  select_population <- pop %>% select(name, "2020")

  # some countries names have names that differ than the ones used in the JHU data
  # This step changes names that differ to the names used in the JHU data
  colnames(select_population)[1] <- "country"
  select_population$country <- as.character(select_population$country)
  select_population$country[str_which(select_population$country, "Korea")][2] <- "South Korea"
  select_population$country[str_which(select_population$country, "Iran")] <- "Iran"
  select_population$country[str_which(select_population$country, "Russia")] <- "Russia"
  select_population$country[str_which(select_population$country, "United States of ")] <- "US"
  select_population$country[str_which(select_population$country, "United Kingdom")] <- "UK"
  select_population$country[str_which(select_population$country, "Taiw")] <- "Taiwan"
  select_population$country[str_which(select_population$country, "Hong")] <- "Hong Kong"



  death_threshold_iter_sample <- mclapply(1:100, function(death_thres) {
    death_threshold_specific_corr(death_thres, 8, countryEMPscore = countryEMPscore_sample)
  }, mc.cores = 10)



  # create data frame from death threshold correlations
  # combine results into a matrix
  df_Death_sample <- do.call(rbind, death_threshold_iter_sample)

  # normalize days for visualization purposes
  norm_df_Death_sample <- do.call(rbind, lapply(death_threshold_iter_sample, function(d) {
    d$days <- ((d$days) - min(d$days)) / ((max(d$days) - min(d$days)))
    d
  }))


  # create matrix of only significant correlations
  # rename corr types
  norm_df_Death_sample$source <- as.character(norm_df_Death_sample$source)
  norm_df_Death_sample <- norm_df_Death_sample[which(norm_df_Death_sample$source == "Structural_proteins"), ]
  norm_df_Death_sample$source[which(norm_df_Death_sample$source == "Structural_proteins")] <- "SARS-CoV-2 structural proteins"
  norm_df_Death_sample$source <- as.factor(norm_df_Death_sample$source)
  sig <- subset(norm_df_Death_sample, sig == 1)

  resample_dist[[rev_2]] <-data.frame(norm_df_Death_sample ,iter=rev_2)
}


KLD_data <- data.frame(SP = sapply(resample_dist, function(i) {

  KLD(px = i$correlation, py = df_death_comp$correlation[which(df_death_comp$source == "Structural_proteins")])$sum.KLD.px.py
}), AP = sapply(resample_dist, function(i) {

  KLD(px = i$correlation, py = df_death_comp$correlation[which(df_death_comp$source == "All_proteins")])$sum.KLD.px.py
})) %>% reshape2::melt()


KLD_box <- KLD_data %>% ggplot(aes(y = value, x = variable,fill=variable)) +
  geom_boxplot(notch = T) +
  coord_cartesian(ylim = c(0, .025)) +
  ylab("KLD") +
  xlab("") +
  theme(legend.position = "bottom")+coord_flip()

wilcox.test(KLD_data$value[which(KLD_data$variable=="AP")],KLD_data$value[which(KLD_data$variable=="SP")])
df <- do.call(rbind, resample_dist)

density_overlay <- df %>% ggplot(aes(x=correlation, group = factor(iter))) +
  geom_density(colour = alpha("grey", .25)) +
  geom_density(inherit.aes = F,data = subset(df_death_comp, source == "Structural_proteins"), aes(correlation), color = "black") +
  geom_density(inherit.aes = F,data = subset(df_death_comp, source == "All_proteins"), aes(correlation), color = "red")+
  geom_density(inherit.aes = F,data = med_cor,aes(med),color="green")


