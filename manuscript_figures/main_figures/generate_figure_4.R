# read in paths
source("~/Covid-19/EnsembleMHC-Covid19/manuscript_figures/set_paths.R")

library(ggpubr)
library(ggplot2)
library(data.table)
library(stringr)
library(wpp2019)
library(dplyr)
library(parallel)
library(lubridate)
library(viridis)
library(pwr)
library(reshape2)


load(paste0(data_path, "selected_countries.R"))
data("pop")

#-------------------------------------------------------
# functions
#-------------------------------------------------------

# min normalize function
min_norm <- function(x) {
  (x - min(x)) / (max(x) - min(x))
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
          })
        })), 2, sum) / length(df$HLA)
      })
    ),
    stringsAsFactors = F
  )
}

# calculate correlation as a function of time
death_threshold_specific_corr_covariate <- function(death_threshold, MIN_countries, method) {
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
    slice(which(Deaths >= death_threshold)) %>%
    merge(corr_data)

  cor_data <- do.call(rbind, lapply(1:thres, function(i) {
    # create subset of data at day i
    tmp <- score_and_death_pop %>% slice(which(days == i))
    dpp <- tmp$death_per_pop
    tmp <- tmp %>% select(-All.proteins, -death_per_pop, -Deaths, -days)
    t(apply(tmp[, -1], 2, function(w) {
      corr_test <- cor.test(dpp, w, method = method)
      c(p_value = corr_test$p.value, corr = corr_test$estimate, day = i, deaths = death_threshold, num_c = length(dpp))
    }))
  })) %>% data.frame()

  cor_data$covar <- row.names(cor_data)[str_which(row.names(cor_data), "[A-Za-z]$")]

  cor_data %>%
    group_by(covar) %>%
    mutate(p_value = p.adjust(p_value, method = "fdr"))
}

# fix the data names
fix_grab_data_names <- function(x, sel) {
  x[str_which(x[, 1], "Iran"), 1] <- "Iran"
  x[str_which(x[, 1], "Russia"), 1] <- "Russia"
  x[str_which(x[, 1], "^Republic of Korea"), 1] <- "South Korea"
  x[str_which(x[, 1], "United King"), 1] <- "UK"
  x[str_which(x[, 1], "United State"), 1] <- "US"
  x[which(x[, 1] %in% sel), ]
}

#-------------------------------------------------
# generate EMP score data
#-------------------------------------------------

# this is taken from the JHU github page: https://github.com/CSSEGISandData/COVID-19
# load JHU data 1/22 - 4/10
all_data_clean <- fread(paste0(data_path, "JHU_data_clean.csv"))
# all_data_clean<-fread("~/Covid-19/Working_dir/improtant_intermediate_datasets/JHU_6.19.csv")


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

# select countries
countryEMPscore <- countryEMPscore[-which(countryEMPscore$country %in% c("Taiwan", "Hong Kong")), ]


#-------------------------------------------------
# covariate data
#-------------------------------------------------

# covariate data

# Prevalence of obesity among adults, BMI ≥ 30, age-standardized
# https://apps.who.int/gho/data/node.main.A900A?lang=en

# Prevalence of overweight among adults, BMI ≥ 25, age-standardized
# https://apps.who.int/gho/data/node.main.A897A?lang=en

# Mean BMI
# same website

# Current expenditure on health by general government and compulsory schemes (% of current expenditure on health)
# https://apps.who.int/gho/data/node.main.HS05?lang=en

# death by non communicable diseases
# https://apps.who.int/gho/data/view.main.NCDDEATHCAUSESNUMBERv?lang=en

# there is no data for Hong knong or for taiwan so those are excluded

#------------------------------
# obesity data age standardized BMI > 30
#------------------------------

obesity <- fread(paste0(data_path, "population_covariates/Obseity_age_standardized.csv")) %>% data.frame()

# fix diff country names
obesity <- fix_grab_data_names(obesity, country_names)

# extract percentage of obese population
obesity$Both.sexes <- as.numeric(str_extract(obesity$Both.sexes, ".+(?= \\[)"))

# rename columns
colnames(obesity) <- c("Country", "percentage_obese")

#------------------------------
# overweight data age standardized BMI >25
#------------------------------
overweight <- fread(paste0(data_path, "population_covariates/Overweight_age_standardized.csv")) %>% data.frame()

# fix diff country names
overweight <- fix_grab_data_names(overweight, country_names)

# extract percentage of overweight population
overweight$Both.sexes <- as.numeric(str_extract(overweight$Both.sexes, ".+(?= \\[)"))

# rename columns
colnames(overweight) <- c("Country", "percentage_overweight")

#------------------------------
# mean BMI
#------------------------------
mean_BMI <- fread(paste0(data_path, "population_covariates/mean_BMI_by_country.csv")) %>% data.frame()

# fix diff country names
mean_BMI <- fix_grab_data_names(mean_BMI, country_names)

# extract mean BMI
mean_BMI$Both.sexes <- as.numeric(str_extract(mean_BMI$Both.sexes, ".+(?= \\[)"))

# rename columns
colnames(mean_BMI) <- c("Country", "avg_BMI")

#--------------------------------------------
# health expenderature as a percentage of GDP
#--------------------------------------------
health_per_GPD <- fread(paste0(data_path, "population_covariates/health_expend_as_percentageGDP.csv")) %>% data.frame()

# fix diff country names
health_per_GPD <- fix_grab_data_names(health_per_GPD, country_names)

# rename columns
colnames(health_per_GPD) <- c("Country", "GDP_expend")

#--------------------------------------------
# age
#--------------------------------------------
age_dist <- rbind(
  popF %>% mutate(age_group = as.numeric(str_extract(popF$age, pattern = "[0-9]+"))),
  popM %>% mutate(age_group = as.numeric(str_extract(popM$age, pattern = "[0-9]+")))
) %>%
  mutate(risk = if_else(age_group >= 65, 1, 0)) %>%
  select(name, pop = `2020`, age_group, risk) %>%
  group_by(name, risk) %>%
  summarise(pop = sum(pop)) %>%
  group_by(name) %>%
  mutate(total_pop = sum(pop)) %>%
  mutate(gt_sixityfive = pop / total_pop) %>%
  select(name, risk, gt_sixityfive)

colnames(age_dist)[1] <- "country"
age_dist$country <- as.character(age_dist$country)
age_dist$country[str_which(age_dist$country, "Korea")][2] <- "South Korea"
age_dist$country[str_which(age_dist$country, "Iran")] <- "Iran"
age_dist$country[str_which(age_dist$country, "Russia")] <- "Russia"
age_dist$country[str_which(age_dist$country, "United States of ")] <- "US"
age_dist$country[str_which(age_dist$country, "United Kingdom")] <- "UK"

age_dist <- age_dist %>%
  ungroup() %>%
  slice(which(country %in% countryEMPscore$country)) %>%
  slice(which(risk == 1)) %>%
  select(-risk, Country = country, gt_sixityfive)

#--------------------------------------------
# sex
#--------------------------------------------
sex_dist <- rbind(
  popF %>% group_by(name) %>% summarise(pop = sum(`2020`)) %>% mutate(sex = "F"),
  popM %>% group_by(name) %>% summarise(pop = sum(`2020`)) %>% mutate(sex = "M")
) %>%
  group_by(name) %>%
  mutate(total = sum(pop)) %>%
  mutate(prop = pop / total) %>%
  select(name, sex, prop) %>%
  slice(which(sex == "M"))

colnames(sex_dist)[1] <- "country"
sex_dist$country <- as.character(sex_dist$country)
sex_dist$country[str_which(sex_dist$country, "Korea")][2] <- "South Korea"
sex_dist$country[str_which(sex_dist$country, "Iran")] <- "Iran"
sex_dist$country[str_which(sex_dist$country, "Russia")] <- "Russia"
sex_dist$country[str_which(sex_dist$country, "United States of ")] <- "US"
sex_dist$country[str_which(sex_dist$country, "United Kingdom")] <- "UK"

sex_dist <- sex_dist %>%
  ungroup() %>%
  slice(which(country %in% countryEMPscore$country)) %>%
  select(Country = country, sex_prop = prop)


#-----------------------------------------------------------------------------------------
# General government expenditure on health as a percentage of total government expenditure
#-----------------------------------------------------------------------------------------
health_expend <- fread(paste0(data_path, "population_covariates/expend_by_gen_gov_and_comp_scheme_per.csv")) %>% data.frame()

# fix diff country names
health_expend <- fix_grab_data_names(health_expend, country_names)

# extract percentage of blood glucose population
health_expend$General.government.expenditure.on.health.as.a.percentage.of.total.government.expenditure <- as.numeric(health_expend$General.government.expenditure.on.health.as.a.percentage.of.total.government.expenditure)

# rename columns
colnames(health_expend) <- c("Country", "perc_expend")

#--------------------------------------------
# deaths by different metrics
#--------------------------------------------
ncom <- fread(paste0(data_path, "population_covariates/death_by_noncomm_dis.csv")) %>%
  data.frame() %>%
  slice(which(Year == 2016))

# fix diff country names
ncom <- fix_grab_data_names(ncom, country_names)

ncom_lst <- lapply(unique(ncom$Causes), function(w) {
  tmp <- ncom[which(ncom$Causes == w), ]
  tmp <- tmp[, c(1, 4)]
  colnames(tmp) <- c("Country", str_replace_all(w, pattern = " ", "_"))
  tmp
})


ncom_mat <- Reduce(function(...) merge(..., by = c("Country")), ncom_lst)

# load population data
data("pop")

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

# process data
# convert to deaths per million due to diff diseases
# get rid of malignant neoplasms because we dont care about these
ncom_proc <- merge(ncom_mat, select_population %>% select(Country = country, pop = `2020`)) %>% mutate(pop = pop / 1000)
non_com_data <- ncom_proc %>%
  mutate(
    Malignant_neoplasms = Malignant_neoplasms / pop,
    Diabetes_mellitus = Diabetes_mellitus / pop,
    Cardiovascular_diseases = Cardiovascular_diseases / pop,
    Chronic_obstructive_pulmonary_disease = Chronic_obstructive_pulmonary_disease / pop
  ) %>%
  select(-pop, -Malignant_neoplasms)
# select(-pop)

#-----------------------------------------------------------------------------------------
# Blood pressure
#-----------------------------------------------------------------------------------------
Blood_pressure <- fread(paste0(data_path, "population_covariates/BP_by_country.csv")) %>% data.frame()

# fix diff country names
Blood_pressure <- fix_grab_data_names(Blood_pressure, country_names)

# extract percentage of blood pressure
Blood_pressure$Raised.blood.pressure..SBP.gt..140.OR.DBP.gt..90...age.standardized.estimate. <- as.numeric(str_extract(Blood_pressure$Raised.blood.pressure..SBP.gt..140.OR.DBP.gt..90...age.standardized.estimate., ".+(?= \\[)"))


# rename columns
colnames(Blood_pressure) <- c("Country", "blood_pressure")

#------------------------------
# combine the datasets into one
#------------------------------

# create the list of data sets
list_of_data <- list(non_com_data, health_per_GPD, mean_BMI, overweight, obesity, health_expend, age_dist, Blood_pressure, sex_dist)

# merge data sets
datasets <- Reduce(function(...) merge(..., by = c("Country")), list_of_data)

# combined correlation data
corr_data <- merge(datasets, countryEMPscore %>% select(Country = country, Structural.proteins))


colnames(corr_data)[1] <- "country"

#--------------------------------------------------------------------
# correlations with individual metrics
#--------------------------------------------------------------------

# find correlations
death_threshold_iter_cors <- mclapply(1:100, function(death_thres) {
  death_threshold_specific_corr_covariate(death_thres, 8, "spearman")
}, mc.cores = 10)


# normalize days cors
df_death_cor <- do.call(rbind, lapply(death_threshold_iter_cors, function(i) {
  i$day <- min_norm(i$day)
  i
})) %>% data.frame()

# label significant entries
df_death_cor$sig <- if_else(df_death_cor$p_value <= .05, 1, 0)

# factor the covariates
df_death_cor$covar <- factor(df_death_cor$covar)

# graph the correlates
cors <- df_death_cor %>% ggplot(aes(x = day, corr.rho, color = deaths, group = deaths)) +
  geom_line() +
  scale_color_viridis() +
  facet_wrap(covar ~ ., nrow = 3) +
  geom_point(data = subset(df_death_cor, sig == 1), aes(x = day, y = corr.rho), color = "red") +
  theme_pubclean()

cors

tab <- table(df_death_cor$covar, df_death_cor$sig)[, 2] / table(df_death_cor$covar)

# best combination models -------------------------------------------------

top_models <- fread(paste0(data_path, "population_covariates/best_model_data.csv"))


df <- top_models %>%
  group_by(ident) %>%
  summarise(med = median(r2), sig = table(p_value <= .05)[2] / sum(table(p_value <= .05))) %>%
  arrange(med) %>%
  melt(id.vars = c("ident"))

df$ident <- factor(df$ident, levels = df$ident[which(df$variable == "med")][order(df$value[which(df$variable == "med")], decreasing = T)])



p2<-df %>%
  ggplot(aes(x = ident, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme(axis.text.x = element_text(angle = 90))

p2