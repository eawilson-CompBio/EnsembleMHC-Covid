a<-commandArgs(trailingOnly = T)
library(data.table)
library(dplyr)

tab<-lapply(list.files(pattern="out"),function(w){
    fread(w,header=T)%>%unique()
})

combo<-Reduce(merge,tab)
combo<-combo%>%unique()

write.csv(combo,file=paste(a,"scored_peptides.csv",sep="_"),row.names=F)

