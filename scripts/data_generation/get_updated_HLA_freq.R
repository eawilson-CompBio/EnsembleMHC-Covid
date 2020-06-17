#----------------------------
# this will get individual HLA data from allelefrequncy.net 
#This very good allele frequncy net database scrapper was traken from https://github.com/pdxgx/covid19/blob/master/src/HLA_frequencies.R
#lines 18-49 are a direct copy from this script 
#However a conversion error in this scrappr was identified where alleles with very low expression <=.0008 were set to errounous values (e.g. 804,404 etc)
#As it would be impossible to achieve allele frequncies greater than 1, all values greater than one are assumed to be casued by this error and are set to 0
#this will not casue a noticible change in EnsembleMHC score calculation as theses values are very small.
#----------------------------

##set path to ensembleMHC-Covid directory
##make sure it does not end with a "/" (i.e. ~/EnsembleMHC-Covid19 )
Ensemble_PATH="~/Covid-19/EnsembleMHC-Covid19"
dataset_path<-paste0(Ensemble_PATH,"/datasets")

#liad libraries
library(rvest)
library(stringr)
library(dplyr)



alleles <- data.frame(HLA=character(), pop.ID=numeric(), pop.name=character(), freq=numeric(), sample.size=numeric())
for (HLA in LETTERS[1:3]) {
  page <- 1
  while (page > 0) {
    print(paste0("HLA-",HLA," (page ",page,")"))
    qurl <- paste0("http://www.allelefrequencies.net/hla6006a.asp?", "hla_locus=", HLA, "&page=", page)
    data <- qurl %>%
      read_html()
    recs <- data %>%
      html_nodes(xpath = '//*[@id="divGenNavig2"]/table') %>%
      html_table()
    recs <- recs[[1]][1]
    if (sub(".*to ([0-9,]+)[^0-9]*[(]from (.*)[)].*", "\\1", recs) ==
        sub(".*to ([0-9,]+)[^0-9]*[(]from (.*)[)].*", "\\2", recs)) {
      page <- 0
    } else {
      page <- page + 1
    }
    data <- data %>%
      html_nodes(xpath = '//*[@id="divGenDetail"]/table') 
    popIDs <- data %>% 
      html_nodes("a") %>% 
      html_attr("href") %>% 
      grep(pattern="pop6001c", value=TRUE) %>%
      sub(pattern=".*[=]",replacement="") %>%
      as.numeric()
    data <- data %>% 
      html_table() %>%
      unlist(recursive=FALSE)
    alleles <- rbind(alleles, data.frame(HLA=data[[2]], pop.ID=popIDs, pop.name=data[[4]], freq=sub("[^0-9,.]+","",data[[6]]), sample.size=gsub(",","",data[[8]])))
  }
}

#grab at least 4 digit HLA
alleles<-alleles[str_which(alleles$HLA,"[A-C]\\*[0-9]+\\:[0-9]+"),]
#extract the 4 digit HLA
alleles$HLA<-str_extract(alleles$HLA,"[A-C]\\*[0-9]+\\:[0-9]+")
#remove *
alleles$HLA<-str_remove(alleles$HLA,"\\*")
#remove comma from sample sizes >999 and convert to numeric
alleles$sample.size<-as.numeric(str_remove(alleles$sample.size,","))
#convert freq to numeric as it is  converted to a factor at some point above
alleles$freq<-as.numeric(levels(alleles$freq)[as.numeric(alleles$freq)])
#set problematic frequncies to 0
alleles$freq[which(alleles$freq>1)]<-0
#select coloumns of interest
y<-alleles%>%select(HLA,pop.name,freq,sample.size)
#rename columns
colnames(y)<-c("HLA","population","Freq","sample_size")

#save file to dataset folder
write.csv(y,file = paste0(dataset_path,"/HLA_data_updated_",str_extract(Sys.time(),"[0-9]+-[0-9]+-[0-9]+" ),".csv"),row.names = F)
