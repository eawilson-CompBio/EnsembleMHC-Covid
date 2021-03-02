
library(rvest)
library(ggpubr)
library(dplyr)
library(stringr)
library(lubridate)
library(patchwork)
url <- "https://www.worldometers.info/coronavirus/"

todays_data <- html_table(read_html(url))[[1]][, c("Country,Other", "Deaths/1M pop")]
countries=c("China","Hong Kong","Japan","France","Taiwan","Iran","S. Korea","Italy","USA","UK","Netherlands","Germany","Morocco","India","Ireland","Poland","Turkey","Mexico","Russia","Israel","Croatia","Czechia","Romania")
todays_data<-todays_data%>%slice(which(`Country,Other`%in%countries))

todays_data$`Country,Other`[which(todays_data$`Country,Other`=="S. Korea")]<-"South Korea"
todays_data$`Country,Other`[which(todays_data$`Country,Other`=="USA")]<-"US"
colnames(todays_data)<-c("country","deaths_per_million")
todays_data$deaths_per_million<-as.numeric(str_remove_all(todays_data$deaths_per_million,","))

countryEMPscore <- read.csv("populationEMPscore.csv")
data<-todays_data%>%merge(countryEMPscore)

data$struct_rank<-rank(data$Structural.proteins)
data$all_rank<-rank(data$All.proteins)
data$death_rank<-rank(data$deaths_per_million)

data$group<-"no"
data$group[which(data$Structural.proteins>median(data$Structural.proteins))]<-"yes"

cor<-cor.test(data$deaths_per_million,data$Structural.proteins,method = "spearman")
cor_text<-paste0("rho = ",round(cor$estimate,3)," ; p-value = ",formatC(cor$p.value, format = "e", digits = 2))

the_date <- today()
p1 <- data %>% ggplot(aes(x = death_rank, y = struct_rank, label = country)) +
  geom_smooth(formula = y ~ x, method = "lm", color = "black") +
  geom_point(aes(color = group), size = 2) +
  ggrepel::geom_text_repel() +
  theme_classic() +
  geom_hline(yintercept = median(data$struct_rank), linetype = "dashed") +
  xlab("deaths per million rank") +
  ylab("EnsembleMHC population rank") +
  theme(legend.position = "none") +
  ggtitle(paste("correlation as of ",the_date),
          subtitle = cor_text)

mann_text<-paste0("Mann-Whitney U test p-value = ",formatC(wilcox.test(data$deaths_per_million[which(data$group=="yes")],data$deaths_per_million[which(data$group=="no")])$p.value,format = "e",digits = 2))

p2 <- data %>% ggplot(aes(x = group, y = deaths_per_million,fill=group,label=country)) +
  geom_boxplot()+
  geom_text_repel() +
  ylab("deaths per 1M") +
  xlab("EnsembleMHC population rank > median?") +
  theme_pubclean() +
  theme(legend.position = "none") +
  ggtitle("High/Low EMP score stratifictation",subtitle = mann_text)

plot<-p1+p2 

ggsave(plot = plot,filename = paste0("todays_data/",the_date,".pdf"),width = 9,height = 4.5)

