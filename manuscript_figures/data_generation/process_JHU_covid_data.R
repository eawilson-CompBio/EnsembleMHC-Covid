data_path <- "~/Covid-19/EnsembleMHC-Covid19/datasets/"
library(dplyr)
library(tidyr)
library(stringr)

# generate coronavrius case-death data
# this is taken from the JHU github page: https://github.com/CSSEGISandData/COVID-19
# the data is up april 10th
# indicate path to JHU data folder. The below path is an example.
# when using newer data, you will want to manually inspect. The naming scheme in the reported data is not static
case_data <- list.files(path = "~/Covid-19/JHU_data/4.10.corona/COVID-19/csse_covid_19_data/csse_covid_19_daily_reports/", pattern = "csv", full.names = T)
# Just as a heads up, they switched the date notation scheme about 60 days into there tracking. Here I process each indepedently then combine
all_data <- do.call(rbind, lapply(case_data[1:60], function(i) {
  x <- read.csv(i, stringsAsFactors = F)
  tmp <- x[, c(2, 4:6)]
  tmp <- tmp %>%
    group_by(Country.Region) %>%
    replace_na(replace = list(Confirmed = 0, Deaths = 0, Recovered = 0)) %>%
    summarise_each(sum)
  tmp$date <- str_extract(i, "[0-9]+-[0-9]+-2020")
  tmp
}))
# this for the second bit of data
all_data_cont <- do.call(rbind, lapply(case_data[61:length(case_data)], function(i) {
  x <- read.csv(i, stringsAsFactors = F)
  tmp <- x[, c(4, 8:10)]
  tmp <- tmp %>%
    group_by(Country_Region) %>%
    replace_na(replace = list(Confirmed = 0, Deaths = 0, Recovered = 0)) %>%
    summarise_each(sum)
  tmp$date <- str_extract(i, "[0-9]+-[0-9]+-2020")
  tmp
}))
# make sure column names are the same
colnames(all_data_cont) <- colnames(all_data)

# bind the list together
all_d <- data.frame(rbind(all_data, all_data_cont), stringsAsFactors = F)
# convert the data to character
all_d$date <- as.character(all_d$date)

# fix troublesome country names so that they align with the names in the HLA data
all_d$Country.Region[str_which(all_d$Country.Region, "Bahamas")] <- "Bahamas"
all_d$Country.Region[str_which(all_d$Country.Region, "Korea")] <- "South Korea"
all_d$Country.Region[str_which(all_d$Country.Region, "Iran")] <- "Iran"
all_d$Country.Region[str_which(all_d$Country.Region, "China")] <- "China"
all_d$Country.Region[str_which(all_d$Country.Region, "United Kingdom")] <- "UK"
all_d$Country.Region[str_which(all_d$Country.Region, "Taiw")] <- "Taiwan"
all_d$Country.Region[str_which(all_d$Country.Region, "Cz")] <- "Czechia"

# this will address any problems arising from multiple entries from the same day
alldd <- all_d %>%
  group_by(Country.Region, date) %>%
  summarise_each(sum) %>%
  as.data.frame()

# arrange by date
all_data_clean <- alldd %>%
  group_by(Country.Region) %>%
  arrange(date)

write.csv(all_data_clean, file = paste0(data_path, "JHU_data_clean.csv"), row.names = F)