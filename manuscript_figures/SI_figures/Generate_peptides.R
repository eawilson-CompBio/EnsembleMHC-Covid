REVISION_DIRECTORY <- "~/Covid-19/EnsembleMHC-Covid19/revision_requests/"
setwd(paste0(REVISION_DIRECTORY,"Viral_peptide_benchmarking/ref/"))
library(stringr)
library(Biostrings)
library(data.table)
library(dplyr)
#this is HEP-C genome polyprot ; dengue poly prot and gag-pol
viral_seqs <- readAAStringSet(paste0(REVISION_DIRECTORY,"revision_datasets/viral_proteins.fa"))
load(paste0(REVISION_DIRECTORY,"revision_datasets/52_HLA_list.R"))

viral_seqs <- lapply(viral_seqs,function(i){
  str_split(i,pattern = "",simplify = T)
})


names(viral_seqs) <- str_extract(names(viral_seqs),"(?<=\\|).{0,10}(?= )")

peptides<- do.call(rbind, lapply(1:length(viral_seqs), function(i) {
  virus <- viral_seqs[[i]]
  peptides <- do.call("c", mclapply(8:14, function(j) {
    do.call("c", lapply(1:(length(virus) - (j - 1)), function(k) {
      paste(virus[k:(k + (j - 1))], collapse = "")
    }))
  }, mc.cores = 7))
  data.frame(peptide = peptides, gene = names(viral_seqs)[i], length = nchar(peptides))
}))

write.csv(peptides,file = paste0(REVISION_DIRECTORY,"revision_datasets/all_viral_peptides.csv"))
# write the files 
lapply(paste0("HLA-",sel_HLA,"_ready_for_predictions.csv"),function(i){
  write.csv(peptides,file=i,row.names = F)
})

