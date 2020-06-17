file<-commandArgs(trailingOnly = T)
library(data.table)
library(dplyr)
library(stringr)
library(Biostrings)
slice<-dplyr::slice

prot<-readAAStringSet("~/last_chance_covid/PA_ucsc_proteome_259contams_viruses_26tumorShared.fasta")

i=file

tmp<-fread(i)

name <- str_extract(i,"(C|G|M|O).+(?=\\.P)")
print(name)

tmp<-tmp%>%select(sequence,deltaForwardReverseScore,accession_number)%>%mutate(peptide=toupper(sequence))%>%
    arrange(desc(deltaForwardReverseScore))%>%dplyr::slice(which(deltaForwardReverseScore>0))%>%
    dplyr::slice(-which(duplicated(peptide)))%>%dplyr::slice(which(nchar(peptide)>=8&nchar(peptide)<=14))%>%select(peptide,accession_number)


##only look at peptide wiht proteins that map 
  tmp<-tmp[which(tmp$accession_number%in%names(prot)),]


 
 decoy <- do.call("c",mclapply(1:nrow(tmp),function(j){
      j<-unlist(c(tmp[j,]))
      pep_len<-nchar(as.character(j[1]))
      pep_prot<-prot[[j[2]]]
   if(nchar(pep_prot)>113){
   random_indices<-sample(1:(nchar(pep_prot)-(pep_len-1)),100)
   }else{
     random_indices<-1:(nchar(pep_prot)-(pep_len-1))
   }
     prot_chop<-str_split(pep_prot,"",simplify = T)
   decoy_peps<-sapply(random_indices,function(l){
     paste(prot_chop[l:(l+(pep_len-1))],collapse = "")
   })
   
 },mc.cores = 10))
 
combo<-rbind(data.frame(peptide=tmp$peptide,ident="target",stringsAsFactors = F),
      data.frame(peptide=decoy,ident="decoy",stringsAsFactors = F))
 
combo$pep_len<-nchar(combo$peptide)

combo<-combo%>%slice(-which(duplicated(combo$peptide)))%>%mutate(HLA="REPLACE")

write.csv(combo,paste0(name,"_target_decoy.csv"),row.names = F)




