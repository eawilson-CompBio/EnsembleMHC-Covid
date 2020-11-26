# read in paths
source("~/Covid-19/EnsembleMHC-Covid19/manuscript_figures/set_paths.R")
library(ggplot2)
library(data.table)
library(stringr)
library(latex2exp)
library(lubridate)
library(viridis)
library(pwr)
library(dplyr)
library(ggpubr)
library(parallel)


# functions ---------------------------------------------------------------

# create score function that will consider individual algorithms
EMPscore_inv <- function(population, count_mat) {
  data.frame(do.call(rbind, lapply(population, function(df) {
    apply(t(sapply(as.character(str_remove(df$HLA, "\\*")), function(i) {
      sapply(c("All proteins", "Structural proteins"), function(j) {
        # pf is the peptide fraction; count matrix
        pf <- count_mat$pep_frac[which(count_mat$HLA == i & count_mat$source == j)]
        score <- (df$pop_freq[which(df$HLA == i)] * pf)
      })
    })), 2, sum) / length(df$HLA)
  })), stringsAsFactors = F)
}

# min function
min_norm <- function(x) {
  (max(x) - x) / (max(x) - min(x))
}


algos <- c("mhcflurry_affinity_percentile", "mhcflurry_presentation_score", "MixMHCpred", "netMHC_affinity", "netMHCpan_EL_affinity", "netstab_affinity", "pickpocket_affinity")

all_data_clean <- fread(paste0(data_path,"/JHU_data_clean.csv"))

corona_data <- fread(paste0(data_path, "predicted_corona_peptides.csv"))

selected_allele <- load(paste0(data_path, "selected_alleles.R"))

# day correction for when the data ends
day_correction <- as.numeric(today() - mdy("4/10/2020"))



load(paste0(data_path, "P_sum_median_1000_boot.R"))
# set arbitrarily low scores
P_sum$value[which(P_sum$algo == "netstab_affinity")] <- .5
P_sum$value[which(P_sum$algo == "netMHCpan_EL_affinity")] <- .5
P_sum$value[which(P_sum$algo == "MixMHCpred")] <- .5
P_sum$value[which(P_sum$algo == "mhcflurry_affinity_percentile")] <- .5
P_sum$value[which(P_sum$algo == "netMHC_affinity")] <- .5
# pickpocket set to 50nm threshold
P_sum$value[which(P_sum$algo == "pickpocket_affinity")] <- 1 - log(50, base = 50000)
# presentation set to the top 0.5 percent
P_sum$value[which(P_sum$algo == "mhcflurry_presentation_score")] <- quantile(corona_data$mhcflurry_presentation_score,.995)


pep_probs_ind <- lapply(unique(corona_data$HLA), function(w) {
  tmp <- corona_data %>% dplyr::slice(which(corona_data$HLA == w)) %>% data.frame()
  tmp$mhcflurry_presentation_score <- min_norm(tmp$mhcflurry_presentation_score)
  tmp$gene <- factor(tmp$gene)
  print(w)
  prob_list <- lapply(colnames(tmp)[colnames(tmp) %in% unique(P_sum$algo)], function(q) {
    if (q == "pickpocket_affinity") {
      neg <- 1 - P_sum$PPV[which(P_sum$HLA == str_remove(w, pattern = "\\*") & P_sum$algo == q)]
      thres <- P_sum$value[which(P_sum$HLA == str_remove(w, pattern = "\\*") & P_sum$algo == q)]
      peptides <- tmp$peptide[which(tmp[, q] >= thres)]
      if (length(peptides) > 1) {
        tmp <- data.frame(peptide = peptides, prob = neg, algo = q)
      }
    } else {
      neg <- 1 - P_sum$PPV[which(P_sum$HLA == str_remove(w, pattern = "\\*") & P_sum$algo == q)]
      thres <- P_sum$value[which(P_sum$HLA == str_remove(w, pattern = "\\*") & P_sum$algo == q)]
      if (is.na(neg)) {
        print(q)
      }
      peptides <- tmp$peptide[which(tmp[, q] <= thres)]
      if (length(peptides) > 1) {
        tmp <- data.frame(peptide = peptides, prob = neg, algo = q)
      }
    }
  })
  prob_combo <- do.call(rbind, prob_list) %>%
    data.frame() %>%
    group_by(peptide) %>%
    merge(tmp[, c("peptide", "gene", "HLA")])
})



all_ind <- do.call(rbind, pep_probs_ind)
all_ind <- all_ind %>% mutate(PPV = 1 - prob)



indv_metric <- lapply(algos, function(w) {
  struct <- all_ind[which(all_ind$algo == w), ] %>%
    slice(which(gene %in% c("E", "N", "M", "S"))) %>%
    group_by(HLA) %>%
    summarise(count = length(peptide)) %>%
    mutate(pep_frac = count / sum(count)) %>%
    mutate(source = "Structural proteins")
  all <- all_ind[which(all_ind$algo == w), ] %>%
    group_by(HLA) %>%
    summarise(count = length(peptide)) %>%
    mutate(pep_frac = count / sum(count)) %>%
    mutate(source = "All proteins")
  rbind(all, struct)
})

# name the indv_metric list
names(indv_metric) <- algos
# check lengths becasue high thresholds might cause some alleles to have no predicted peptides
algo_lengths <- sapply(indv_metric, nrow)


indv_metric[[which(algo_lengths != 104)]] <- do.call(rbind, lapply(c("Structural proteins", "All proteins"), function(j) {
  tmp <- indv_metric[[which(algo_lengths != 104)]]
  tmp <- tmp[which(tmp$source == j), ]
  rbind(tmp, data.frame(HLA = sel_alleles[-which(sel_alleles %in% tmp$HLA)], count = 0, pep_frac = 0, source = j))
})) 


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




indv_metric <- c(indv_metric, list(HLA_count))
names(indv_metric) <- c(algos, "EnsembleMHC")



all_HLA_data <- fread(paste0(data_path, "HLA_data_manuscript.csv"), header = T)
# countries with at least 1 case reported case by 04/10
test_set <- unique(all_data_clean$Country.Region[which(all_data_clean$Deaths >= 1)])
# save country names
country_list <- test_set

# fix country names so that they will match HLA data
test_set <- lapply(test_set, function(w) {
  w
})
test_set[[which(country_list == "Czechia")]] <- "Czech"
test_set[[which(country_list == "UK")]] <- c("United Kingdom|England|Wales$")
test_set[[which(country_list == "US")]] <- "USA"



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
# get corona virus data that correspondes to the selected alleles
CoV_data_sel_countries <- all_data_clean[which(all_data_clean$Country.Region %in% names(allele_freq_list)), ]


plots <- lapply(1:length(indv_metric), function(n) {
  score_mat <- indv_metric[[n]]
  # apply the peptide score function
  countryEMPscore <- EMPscore_inv(allele_freq_list, score_mat)
  # add a column for the names of the countries
  countryEMPscore$country <- row.names(countryEMPscore)

  # fix teh column name so you can merge later
  colnames(CoV_data_sel_countries)[1] <- "country"


  library(wpp2019)

  data("pop")

  # load pop data and fix names
  select_population <- pop %>% select(name, "2020")

  colnames(select_population)[1] <- "country"
  select_population$country <- as.character(select_population$country)
  select_population$country[str_which(select_population$country, "Korea")][2] <- "South Korea"
  select_population$country[str_which(select_population$country, "Iran")] <- "Iran"
  select_population$country[str_which(select_population$country, "Russia")] <- "Russia"
  select_population$country[str_which(select_population$country, "United States of ")] <- "US"
  select_population$country[str_which(select_population$country, "United Kingdom")] <- "UK"
  select_population$country[str_which(select_population$country, "Taiw")] <- "Taiwan"
  select_population$country[str_which(select_population$country, "Hong")] <- "Hong Kong"

  death_threshold_specific_corr <- function(death_threshold, MIN_countries) {
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


    #  melt_df_pre$confirmed<-death_threshold
    melted_df$confirmed <- death_threshold
    # melt_df_pre
    # return matrix
    melted_df
  }



  death_threshold_iter <- mclapply(1:100, function(death_thres) {
    death_threshold_specific_corr(death_thres, 8)
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
  norm_df_Death$source[which(norm_df_Death$source == "Structural_proteins")] <- "SARS-CoV-2 structural proteins"
  norm_df_Death$source[which(norm_df_Death$source == "All_proteins")] <- "entire SARS-CoV-2 proteome"
  norm_df_Death$source <- as.factor(norm_df_Death$source)
  sig <- subset(norm_df_Death, sig == 1)
  print(table(sig$source) / table(norm_df_Death$source))


  df_spearman_pwr <- norm_df_Death
  df_spearman_pwr$pwr <- apply(df_spearman_pwr, 1, function(w) {
    pwr.r.test(r = as.numeric(w[3]), n = as.numeric(w[6]), sig.level = .05)[[4]]
  })

  R <- 1
  df_spearman_pwr$PPV <- (df_spearman_pwr$pwr * R) / ((df_spearman_pwr$pwr * R) + df_spearman_pwr$p_value)
  ggplot(norm_df_Death, aes(days, correlation, group = confirmed)) +
    geom_line(aes(color = confirmed)) +
    scale_color_viridis() +
    geom_point(data = subset(df_spearman_pwr, PPV >= .95), aes(days, correlation), color = "red", size = .25) +
    facet_wrap(source ~ .) +
    theme_linedraw() +
    theme(legend.position = "bottom") +
    labs(color = "number of deaths by day 0") +
    guides(colour = guide_colorbar(title.position = "bottom")) +
    ylab(TeX("$correlation(\\rho)$")) +
    xlab("normalized days") +
    ggtitle(names(indv_metric)[n])
})


garg <- do.call(ggarrange, plots)

garg
#ggsave(garg, filename = paste0(EnsembleMHC_path, "plots/SI_figures/SI_individual_algorithms_PPV.pdf"), height = 12, width = 15)
