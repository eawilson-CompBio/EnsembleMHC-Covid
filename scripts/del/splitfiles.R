file<-commandArgs(trailingOnly = T)
library(data.table)
library(dplyr)
library(stringr)
library(Biostrings)


prot<-readAAStringSet("PA_ucsc_proteome_259contams_viruses_26tumorShared.fasta")

i=file

tmp<-fread(i)

HLA=paste("HLA-",str_extract(i,"[A-C][0-9]{2}"),":",str_extract(i,"(?<=[A-C][0-9]{2})[0-9]{2}"),sep="")

tmp<-tmp%>%select(sequence,deltaForwardReverseScore,accession_number)%>%mutate(peptide=toupper(sequence))%>%
    arrange(desc(deltaForwardReverseScore))%>%dplyr::slice(which(deltaForwardReverseScore>0))%>%
    dplyr::slice(-which(duplicated(peptide)))%>%dplyr::slice(which(nchar(peptide)>7&nchar(peptide)<16))%>%select(peptide,accession_number)%>%
    mutate(HLA=HLA)


##only look at peptide wiht proteins that map 
  tmp<-tmp[which(tmp$accession_number%in%names(prot)),]

target <- t(apply(tmp,1,function(j){

    pep<-str_extract(prot[[as.character(j[2])]],paste(".{0,3}",j[1],".{0,3}",sep = ""))
    
    if(nchar(pep)<nchar(j[1])+6){
      split<-str_split(pep,pattern=j[1],simplify = T)
      for(k in 1:2){
        if(nchar(split[k])<3&k==1){
          split[1]<-trimws(paste(paste(rep("-",3-nchar(split[1])),collapse = ""),split[1],sep = ""),"both")
        }
        if(nchar(split[k])<3&k==2){
          split[2]<-trimws(paste(paste(rep("-",3-nchar(split[2])),collapse = ""),split[2],sep = ""),"both")
        }
      }
    }else{
      split<-str_split(pep,pattern=j[1],simplify = T)
    }
    
    c(j[[1]],str_split(paste(split,collapse = ""),"",simplify = T),as.character(j[2]),nchar(prot[[as.character(j[2])]]),nchar(j[1]))
    
    
  }))
 
 decoy <-  do.call(rbind,lapply(1:nrow(tmp),function(j){
      j<-unlist(c(tmp[j,]))
      pep_len<-nchar(j[1])
   pep_prot<-prot[[as.character(j[2])]]
   if(nchar(pep_prot)>110){
   random_indices<-sample(1:(nchar(pep_prot)-(pep_len-1)),100)
   }else{
     random_indices<-1:(nchar(pep_prot)-(pep_len-1))
   }
     prot_chop<-str_split(pep_prot,"",simplify = T)
   decoy_peps<-sapply(random_indices,function(l){
     paste(prot_chop[l:(l+(pep_len-1))],collapse = "")
   })
   t(sapply(decoy_peps,function(m){
     
   pep<-str_extract(pep_prot,paste(".{0,3}",m[1],".{0,3}",sep = ""))
   
   if(nchar(pep)<nchar(m[1])+6){
     split<-str_split(pep,pattern=m[1],simplify = T)
     for(k in 1:2){
       if(nchar(split[k])<3&k==1){
         split[1]<-trimws(paste(paste(rep("-",3-nchar(split[1])),collapse = ""),split[1],sep = ""),"both")
       }
       if(nchar(split[k])<3&k==2){
         split[2]<-trimws(paste(paste(rep("-",3-nchar(split[2])),collapse = ""),split[2],sep = ""),"both")
       }
     }
   }else{
     split<-str_split(pep,pattern=m[1],simplify = T)
   }
   
   c(m[[1]],str_split(paste(split,collapse = ""),"",simplify = T),as.character(j[2]),nchar(pep_prot),nchar(j[1]))
   
   }))
   
 }))
 
 colnames(target)<-c("peptide",paste("up_",1:3,sep=""),paste("down_",1:3,sep=""),"ACC","prot_length","pep_len")
 colnames(decoy)<-c("peptide",paste("up_",1:3,sep=""),paste("down_",1:3,sep=""),"ACC","prot_length","pep_len")

 tars<-data.frame(target,ident="target",stringsAsFactors = F)
 decs<-data.frame(decoy,ident="decoy",stringsAsFactors = F)
 peptides<-rbind(tars,decs)

 peptides<-peptides[which(apply(peptides,1,function(w){ sum(w=="X") })==0),]
 
 out<-data.frame(peptides,HLA=HLA,stringsAsFactors = F)
 write.csv(out,paste("./pred_dir/",HLA,"_ready_for_prediction.csv",sep=""),row.names=F)


