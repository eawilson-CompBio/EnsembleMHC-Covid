library(rvest)
library(parallel)
library(stringr)
data_path <- "~/Covid-19/EnsembleMHC-Covid19/datasets/"
# This HLA scrapper was heavily adapted from
# Nguyen, Austin, et al. "Human leukocyte antigen susceptibility map for SARS-CoV-2." Journal of virology (2020).
# (https://github.com/pdxgx/covid19/blob/master/src/HLA_frequencies.R)
# major changes are this runs multicore and is focused on pulling the data needed for EnsembleMHC population score
# if you would like to see the difference, please check out https://github.com/pdxgx/covid19/blob/master/src/HLA_frequencies.R

HLA_list <- list()
length(HLA_list) <- 3
names(HLA_list) <- LETTERS[1:3]
for (HLA in LETTERS[1:3]) {
  # get the total number of pages
  qurl <- paste0("http://www.allelefrequencies.net/hla6006a.asp?", "hla_locus=", HLA, "&page=", 1)
  data <- qurl %>%
    read_html()
  recs <- data %>%
    html_nodes(xpath = '//*[@id="divGenNavig2"]/table') %>%
    html_table()
  recs <- recs[[1]][1]
  last_page <- ceiling(as.numeric(str_remove(sub(".*to ([0-9,]+)[^0-9]*[(]from (.*)[)].*", "\\2", recs), ",")) / 100)

  HLA_list[[HLA]] <- do.call(rbind, mclapply(1:last_page, function(page) {
    print(paste(HLA, page))
    qurl <- paste0("http://www.allelefrequencies.net/hla6006a.asp?", "hla_locus=", HLA, "&page=", page)
    data <- qurl %>%
      read_html()

    data <- data %>%
      html_nodes(xpath = '//*[@id="divGenDetail"]/table')

    data <- data %>%
      html_table() %>%
      unlist(recursive = FALSE)

    HLA_data <- data.frame(HLA = as.character(data$Allele), population = as.character(data[[4]]), Freq = as.numeric(str_remove(data[[6]], " \\(\\*\\)")), sample_size = as.numeric(str_remove_all(data[[8]], ",")), stringsAsFactors = F)
    # alleles <- rbind(alleles, data.frame(HLA=data[[2]], pop.ID=popIDs, pop.name=data[[4]], freq=sub("[^0-9,.]+","",data[[6]]), sample.size=gsub(",","",data[[8]])))
  }, mc.cores = 10))
}

All_hla <- do.call(rbind, HLA_list)
All_hla <- All_hla[str_which(All_hla$HLA, "[A-C]\\*[0-9]+\\:[0-9]+"), ]
All_hla$HLA <- str_extract(All_hla$HLA, "[A-C]\\*[0-9]+\\:[0-9]+")

write.csv(All_hla, file = paste(data_path, "HLA_allele_freq_data_", str_split(Sys.time(), pattern = " ", simplify = T)[1], ".csv", sep = ""), row.names = F)
