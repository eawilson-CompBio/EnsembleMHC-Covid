fasta<-commandArgs(trailingOnly = T)
library(stringr)
library(Biostrings)

fasta<-readAAStringSet(fasta)

peptides<-do.call(rbind,lapply(1:length(fasta),function(i){
    peps<-do.call("c",sapply(8:14,function(j){
    prot <- str_split(as.character(fasta[[i]]),pattern = "",simplify = T)
    end  <- length(prot)-(j-1)
      sapply(1:end,function(k){
        paste(prot[k:(k+(j-1))],collapse = "")
      })
  }))
  data.frame(peptide=peps,protein=names(fasta)[i],length=nchar(peps))
}))


write.csv(peptides,file = "EnsembleMHC_pep_pred.tmp",row.names = F)