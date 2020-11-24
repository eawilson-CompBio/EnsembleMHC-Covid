fasta <- commandArgs(trailingOnly = T)
if (!require("stringr")) install.packages("stringr")
if (!require("seqinr")) install.packages("seqinr")
library(stringr)
library(seqinr)

fasta <- read.fasta(fasta, seqtype = "AA")

peptides <- do.call(rbind, lapply(1:length(fasta), function(i) {
  peps <- do.call("c", sapply(8:14, function(j) {
    prot <- as.character(fasta[[i]])
    end <- length(prot) - (j - 1)
    sapply(1:end, function(k) {
      paste(prot[k:(k + (j - 1))], collapse = "")
    })
  }))
  data.frame(peptide = peps, protein = names(fasta)[i], length = nchar(peps))
}))


write.csv(peptides, file = "EnsembleMHC_pep_pred.tmp", row.names = F)
